subroutine readStructuredCGNS(cgns_file)
  use gridData
  use gridInput
  use communication
  use structuredCGNSGrid

  implicit none
  include 'cgnslib_f.h'

  ! Input Arguments
  character*(*),intent(in) :: cgns_file

  ! Working 
  integer(kind=intType) :: i, j, k, ii, istart, iend, localsize, iProc, iZone, offset
  integer(kind=intType):: ierr, base, nZones, nNodes, dims(9), cg, cursize, coorsize
  integer(kind=intType):: CellDim, PhysDim
  integer(kind=intType) :: nbases, start(3), tmpSym, nSymm
  character*100 fileName
  character*32 :: zoneName, bocoName, famName, connectName, donorName, basename
  integer(kind=intType) :: npts, nbocos, bocotype, pts(3, 2), famID, boco, iB2b, transform(3)
  integer(kind=intType) :: donorRange(3, 2), nB2b
  integer(kind=intType) :: ptset_type, normalIndex(3), NormalListFlag, datatype, ndataset
  real(kind=8)   ::  data_double(6), avgNodes, symmSum(3)
  real(kind=realType), dimension(:, :, :), allocatable :: coorX, coorY, coorZ
#ifdef USE_COMPLEX
  complex(kind=realType), dimension(:, :), allocatable :: allNodes, localNodes
#else
  real(kind=realType), dimension(:, :), allocatable :: allNodes, localNodes
#endif
  integer(kind=intType), dimension(:), allocatable :: wallNodes, localWallNodes
  integer(kind=intType), dimension(:, :), allocatable :: sizes

  integer(kind=intType) :: status(MPI_STATUS_SIZE)

  iSymm = 0
  ! Only do reading on root proc:
  if (myid == 0) then 
     print *, ' -> Reading structured CGNS File: ', cgns_file

     ! Open and get the number of zones:
     call cg_open_f(trim(cgns_file), CG_MODE_READ, cg, ierr)
     if (ierr .eq. CG_ERROR) call cg_error_exit_f

     call cg_nbases_f(cg, nbases, ierr)
     if (ierr .eq. CG_ERROR) call cg_error_exit_f
     
     if (nbases .gt. 1) then
        print *, ' ** Warning: pyWarpUstruct only reads the first base in a cgns file'
     end if

     base = 1_intType

     call cg_base_read_f(cg, base, basename, CellDim, PhysDim, ierr)
     if (ierr .eq. CG_ERROR) call cg_error_exit_f
     if (cellDim .ne. 3 .or. PhysDim .ne. 3) then
        print *, 'The Cells must 3 dimensional'
        stop
     end if

     ! !       *** base attribute:  GOTO base node
     ! call cg_goto_f(cg, base, ierr, 'end')
     ! if (ier .eq. CG_ERROR) call cg_error_exit_f

     call cg_nzones_f(cg, base, nZones, ierr);
     if (ierr .eq. CG_ERROR) call cg_error_exit_f
     print *, '   -> Number of Zones:', nzones
     

     allocate(blocks(nZones))
     ! Now count up the total number of nodes
     nNodes = 0
     allocate(sizes(3, nZones))
     do i=1, nZones
        call cg_zone_read_f(cg, base, i, zoneName, dims, ierr)
        sizes(:, i) = dims(1:3)
        nNodes = nNodes + dims(1)*dims(2)*dims(3)
     end do

     ! Now we know the total number of nodes we can allocate the final
     ! required space and read them in. 
     allocate(allNodes(3, nNodes))
     allocate(wallNodes(nNodes))
     wallNodes = 0

     ! Loop back over the nodes and read
     ii = 0;
     familyList(:) = ''
     nWallFamilies = 0
     offset = 0
     zoneLoop: do iZone=1,nZones
        call cg_zone_read_f(cg, base, iZone, zoneName, dims, ierr)
        if (ierr .eq. CG_ERROR) call cg_error_exit_f
        blocks(iZone)%il = dims(1)
        blocks(iZone)%jl = dims(2)
        blocks(iZone)%kl = dims(3)
        blocks(iZone)%name = zoneName
        curSize = dims(1)*dims(2)*dims(3)
        start = (/1, 1, 1/)
        allocate(coorX(dims(1), dims(2), dims(3)), &
                 coorY(dims(1), dims(2), dims(3)), &
                 coorZ(dims(1), dims(2), dims(3)))
        call cg_coord_read_f(cg, base, iZone, "CoordinateX", RealDouble, start, dims(1:3), coorX, ierr)
        if (ierr .eq. CG_ERROR) call cg_error_exit_f
        call cg_coord_read_f(cg, base, iZone, "CoordinateY", RealDouble, start, dims(1:3), coorY, ierr)
        if (ierr .eq. CG_ERROR) call cg_error_exit_f
        call cg_coord_read_f(cg, base, iZone, "CoordinateZ", RealDouble, start, dims(1:3), coorZ, ierr)
        if (ierr .eq. CG_ERROR) call cg_error_exit_f
        
        ! Now interlace the packing      
        do k=1, dims(3)
           do j=1, dims(2)
              do i=1, dims(1)
                 ii = ii + 1
#ifdef USE_COMPLEX
                 allNodes(1, ii) = cmplx(coorX(i, j, k), 0.0)
                 allNodes(2, ii) = cmplx(coorY(i, j, k), 0.0)
                 allNodes(3, ii) = cmplx(coorZ(i, j, k), 0.0)
#else
                 allNodes(1, ii) = coorX(i, j, k)
                 allNodes(2, ii) = coorY(i, j, k)
                 allNodes(3, ii) = coorZ(i, j, k)
#endif
              end do
           end do
        end do

        call cg_nbocos_f(cg, base, iZone, nbocos, ierr)
        if (ierr .eq. ERROR) call cg_error_exit_f
        allocate(blocks(iZone)%bocos(nbocos))

        bocoLoop: do boco=1, nbocos
           call cg_boco_info_f(cg, base, iZone, boco, boconame, bocotype, &
                ptset_type, npts, NormalIndex, NormalListFlag, datatype, ndataset, ierr)
           if (ierr .eq. CG_ERROR) call cg_error_exit_f
           
           call cg_boco_read_f(cg, base, iZone, boco, pts, data_double, ierr)
           if (ierr .eq. CG_ERROR) call cg_error_exit_f
          
           ! Save the boco info
           blocks(iZone)%bocos(boco)%name = boconame
           blocks(iZone)%bocos(boco)%bocoType = bocoType
           blocks(iZone)%bocos(boco)%ptRange = pts
           blocks(iZone)%bocos(boco)%family = ""

           ! Determine if this is in fact a face-bc
           if ( (pts(1,1) == pts(1,2)) .or. (pts(2,1) == pts(2,2)) .or. (pts(3,1) == pts(3,2))) then 
              call cg_goto_f(cg, base, ierr, "Zone_t", iZone, "ZoneBC_t", 1, "BC_t", boco, "end")
              if (ierr == 0) then ! Node exists
                 ! Read family  name if possible
                 call cg_famname_read_f(famName, ierr)
                 !Only if no error:
                 if (ierr .ne. CG_ERROR) then
                    blocks(iZone)%bocos(boco)%family = famName
                 end if

                 if (BCTypeName(bocoType) == 'BCWallViscous' .or. &
                      BCTypeName(bocoType) == 'BCWallInviscid' .or. &
                      BCTypeName(bocoType) == 'BCWall' .or. &
                      BCTypeName(bocoType) == 'BCWallViscousHeatFlux' .or. &
                      BCTypeName(bocoType) == 'BCWallViscousIsothermal') then
                    call checkInFamilyList(familyList, famName, nwallFamilies, famID)
                    if (famID == 0) then
                       nwallFamilies = nwallFamilies + 1
                       familyList(nwallFamilies) = famName
                    end if

                    ! Flag the wall nodes:
                    do k=pts(3, 1), pts(3, 2)
                       do j=pts(2, 1), pts(2, 2)
                          do i=pts(1, 1), pts(1, 2)
                             wallNodes(offset + (k-1)*dims(1)*dims(2) + (j-1)*dims(1) + i) = 1
                          end do
                       end do
                    end do
                 end if
              end if

              if (BCTypeName(bocoType) == 'BCSymmetryPlane') then 
                 ! Determine what the axis of this symmetry plane
                 ! is. We look at the two end points of the range
                 symmSum = zero
                 if (pts(1,1) == pts(1, 2)) then 
                    symmSum(1) = sum(abs(coorX(pts(1,1), :, :)))
                    symmSum(2) = sum(abs(coorY(pts(1,1), :, :)))
                    symmSum(3) = sum(abs(coorZ(pts(1,1), :, :)))
                    nSymm = dims(2)*dims(3)
                 else if (pts(2,1) == pts(2, 2)) then 
                    symmSum(1) = sum(abs(coorX(:, pts(2,1), :)))
                    symmSum(2) = sum(abs(coorY(:, pts(2,1), :)))
                    symmSum(3) = sum(abs(coorZ(:, pts(2,1), :)))
                    nSymm = dims(1)*dims(3)
                 else if (pts(3,1) == pts(3, 2)) then 
                    symmSum(1) = sum(abs(coorX(:, :, pts(3,1))))
                    symmSum(2) = sum(abs(coorY(:, :, pts(3,1))))
                    symmSum(3) = sum(abs(coorZ(:, :, pts(3,1))))
                    nSymm = dims(1)*dims(2)
                 end if

                 if (symmSum(1) / nSymm < symmTol) then 
                    tmpSym = 1
                 else if (symmSum(2) / nSymm < symmTol) then 
                    tmpSym = 2
                 else if (symmSum(3) / nSymm < symmTol) then
                    tmpSym = 3
                 end if
                 
                 if (iSymm == 0) then ! Currently assigned:
                    iSymm = tmpSym
                 end if
               

                 ! Hard zero the sym plane in case our mesh is bad
                 do k=pts(3, 1), pts(3, 2)
                    do j=pts(2, 1), pts(2, 2)
                       do i=pts(1, 1), pts(1, 2)
                          allNodes(isymm, offset + (k-1)*dims(1)*dims(2) + (j-1)*dims(1) + i) = zero
                       end do
                    end do
                 end do

                 if (tmpSym /= iSymm) then
                    print *, 'Error: detected more than 1 symmetry plane direction.'
                    print *, 'pyWarpUnstruct cannot handle this'
                    stop
                 end if
              end if
           end if
        end do bocoLoop
        offset = offset + dims(1)*dims(2)*dims(3)

        call cg_n1to1_f(cg, base, iZone, nB2B, ierr)
        if (ierr .eq. CG_ERROR) call cg_error_exit_f
        allocate(blocks(iZone)%B2Bs(nB2B))

        B2BLoop: do iB2B=1,nB2B
           call cg_1to1_read_f(cg, base, iZone, iB2B, connectName, donorName, pts, &
                donorRange, transform, ierr)
           if (ierr .eq. CG_ERROR) call cg_error_exit_f

           ! Save the b2b info
           blocks(iZone)%b2bs(ib2b)%name = connectName
           blocks(iZone)%b2bs(ib2b)%donorName = donorName
           blocks(iZone)%b2bs(ib2b)%ptRange = pts
           blocks(iZone)%b2bs(ib2b)%donorRange = donorRange
           blocks(iZone)%b2bs(ib2b)%transform = transform
        end do B2BLoop

        ! And free the temporary arrays
        deallocate(coorX, coorY, coorZ)

     end do zoneLoop

     call cg_close_f(cg, ierr)
     if (ierr .eq. CG_ERROR) call cg_error_exit_f
     deallocate(sizes)
  end if

  ! Communication the symmetry direction to everyone
  call MPI_bcast(iSymm, 1, MPI_INTEGER, 0, warp_comm_world, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  
  ! Communicate the total number of nodes to everyone
  call MPI_bcast(nNodes, 1, MPI_INTEGER, 0, warp_comm_world, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  if (myid == 0) then

     istart = 0 + 1
     iend   = int(nNodes*(one/nProc))
     localSize = iend - istart + 1
     allocate(localNodes(3, localSize))
     allocate(localWallNodes(localSize))

     ! Just copy for the root proc:
     localNodes(:, :) = allNodes(:, 1:localSize)
     localWallNodes(:) = wallNodes(1:localSize)

     ! Loop over the remainer of the procs and send
     do iProc=1, nProc-1
        istart = int(nNodes*(dble(iProc)/nProc)) +1 
        iend   = int(nNodes*(dble(iProc+1)/nProc))
        localSize = iend - istart + 1
#ifdef USE_COMPLEX
        call MPI_Send(allNodes(:, iStart), localSize*3, MPI_COMPLEX16, iProc, &
             11, warp_comm_world, ierr)
#else
        call MPI_Send(allNodes(:, iStart), localSize*3, MPI_REAL8, iProc, &
             11, warp_comm_world, ierr)
#endif
        call EChk(ierr, __FILE__, __LINE__)

        call MPI_Send(wallNodes(iStart), localSize, MPI_INTEGER4, iProc, &
             12, warp_comm_world, ierr)
        call EChk(ierr, __FILE__, __LINE__)
     end do
     deallocate(allNodes)
     deallocate(wallNodes)
  else

     istart = int(nNodes*(dble(myid)/nProc)) +1 
     iend   = int(nNodes*(dble(myid+1)/nProc))
     localSize = iend - istart + 1
     
     allocate(localNodes(3, localSize))
     allocate(localWallNodes(localSize))

     ! Receive on all the other procs:
#ifdef USE_COMPLEX
     call MPI_recv(localNodes, 3*localSize, MPI_COMPLEX16, 0, 11, &
          warp_comm_world, status, ierr)
#else
     call MPI_recv(localNodes, 3*localSize, MPI_REAL8, 0, 11, &
          warp_comm_world, status, ierr)
#endif
     call EChk(ierr, __FILE__, __LINE__)
     call MPI_recv(localWallNodes, localSize, MPI_INTEGER4, 0, 12, &
          warp_comm_world, status, ierr)
     call EChk(ierr, __FILE__, __LINE__)
  end if

  ! All we are doing to do here is to create the Xv vector --- which
  ! is done via a call to createCommonGrid

  call createCommonGrid(localNodes, localWallNodes, size(localNodes, 2))
  deallocate(localNodes, localWallNodes)
end subroutine readStructuredCGNS

subroutine checkInFamilyList(familyList, famName, nFamily, index)
  use precision 
  use constants
  implicit none
  character*32, dimension(maxFamilies) :: familyList
  character*32 :: famName
  integer(kind=intType) :: nFamily, i, index

  index = 0_intType
  famLoop: do i=1, nFamily
     if (famName == familyList(i)) then
        index = i
        exit famLoop
     end if
  end do famLoop

end subroutine checkInFamilyList
