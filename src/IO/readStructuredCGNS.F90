subroutine readStructuredCGNS(cg)
  use gridData
  use gridInput
  use communication
  use CGNSGrid

  implicit none
  include 'cgnslib_f.h'

  ! Input Arguments
  integer(kind=intType) :: cg

  ! Working 
  integer(kind=intType) :: i, j, k, ii, istart, iend, localsize, iProc, iZone, offset
  integer(kind=intType):: ierr, base, nZones, nNodes, dims(9),  coorsize
  integer(kind=intType) :: nbases, start(3), tmpSym, nSymm
  character*100 fileName
  character*32 :: zoneName, bocoName, famName
  integer(kind=intType) :: npts, nbocos, bocotype, pts(3, 2), famID, boco
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
  integer(kind=intType) :: nWall, nodeCount, nConn, ni, nj
  integer(kind=intType) :: status(MPI_STATUS_SIZE)
  logical :: lowFace

  iSymm = 0
  ! Only do reading on root proc:
  if (myid == 0) then 

     ! Always the first base
     base = 1_intType
   
     call cg_nzones_f(cg, base, nZones, ierr);
     if (ierr .eq. CG_ERROR) call cg_error_exit_f
     allocate(zones(nZones))
     
     ! Now count up the total number of nodes
     nNodes = 0
     allocate(sizes(3, nZones))
     do i=1, nZones
        call cg_zone_read_f(cg, base, i, zoneName, dims, ierr)
        sizes(:, i) = dims(1:3)
        nNodes = nNodes + dims(1)*dims(2)*dims(3)

        ! Nullify section pointers since we won't have any for structured mesh
        nullify(zones(i)%sections)
     end do

     ! Now we know the total number of nodes we can allocate the final
     ! required space and read them in. 
     allocate(allNodes(3, nNodes))
     allocate(wallNodes(nNodes))
     wallNodes = 0

     ! We make two loops over the zones....the first *just* reads the
     ! BC info and stores it. The second, actually reads the
     ! nodes. 

     familyList(:) = ''
     nWallFamilies = 0
     nWall = 0
     nConn = 0
     zoneLoop_one: do iZone=1,nZones
        call cg_zone_read_f(cg, base, iZone, zoneName, dims, ierr)
        if (ierr .eq. CG_ERROR) call cg_error_exit_f
        zones(iZone)%il = dims(1)
        zones(iZone)%jl = dims(2)
        zones(iZone)%kl = dims(3)
        zones(iZone)%name = zoneName

        call cg_nbocos_f(cg, base, iZone, nbocos, ierr)
        if (ierr .eq. ERROR) call cg_error_exit_f
        allocate(zones(iZone)%bocos(nbocos))

        bocoLoop_one: do boco=1, nbocos

           ! Nullify elemPtr, elemConn and elemNodes since we won't have any
           nullify(zones(iZone)%bocos(boco)%elemPtr, &
                   zones(iZone)%bocos(boco)%elemConn, &
                   zones(iZone)%bocos(boco)%elemNodes) 

           call cg_boco_info_f(cg, base, iZone, boco, boconame, bocotype, &
                ptset_type, npts, NormalIndex, NormalListFlag, datatype, &
                ndataset, ierr)
           if (ierr .eq. CG_ERROR) call cg_error_exit_f
           
           call cg_boco_read_f(cg, base, iZone, boco, pts, data_double, ierr)
           if (ierr .eq. CG_ERROR) call cg_error_exit_f
          
           ! Save the boco info
           zones(iZone)%bocos(boco)%name = boconame
           zones(iZone)%bocos(boco)%type = bocoType
           zones(iZone)%bocos(boco)%ptRange = pts
           zones(iZone)%bocos(boco)%family = ""

           ! Determine if this is in fact a face-bc
           if ( (pts(1,1) == pts(1,2)) .or. (pts(2,1) == pts(2,2)) .or. (pts(3,1) == pts(3,2))) then 
              call cg_goto_f(cg, base, ierr, "Zone_t", iZone, "ZoneBC_t", 1, "BC_t", boco, "end")
              if (ierr == 0) then ! Node exists
                 ! Read family  name if possible
                 call cg_famname_read_f(famName, ierr)
                 !Only if no error:
                 if (ierr .ne. CG_ERROR) then
                    zones(iZone)%bocos(boco)%family = famName
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

                    if (pts(1,1) == pts(1, 2)) then 
                       nWall = nWall + dims(2)*dims(3)
                       nConn = nConn + (dims(2)-1)*(dims(3)-1)
                    else if (pts(2,1) == pts(2, 2)) then 
                       nWall = nWall + dims(1)*dims(3)
                       nConn = nConn + (dims(1)-1)*(dims(3)-1)
                    else if (pts(3,1) == pts(3, 2)) then 
                       nWall = nWall + dims(1)*dims(2)
                       nConn = nConn + (dims(1)-1)*(dims(2)-1)
                    end if
                 end if
              end if

              if (BCTypeName(bocoType) == 'BCSymmetryPlane') then 
                 ! For efficiency, we only load in the points on the 
                 ! symm plane. 

                 if (pts(1,1) == pts(1, 2)) then 
                    allocate(&
                         coorX(pts(2,2)-pts(2,1)+1, pts(3,2)-pts(3,1)+1, 1), &
                         coorY(pts(2,2)-pts(2,1)+1, pts(3,2)-pts(3,1)+1, 1), &
                         coorZ(pts(2,2)-pts(2,1)+1, pts(3,2)-pts(3,1)+1, 1))
                    nSymm = dims(2)*dims(3)
                 else if(pts(2,1) == pts(2,2)) then
                    allocate(&
                         coorX(pts(1,2)-pts(1,1)+1, pts(3,2)-pts(3,1)+1, 1), &
                         coorY(pts(1,2)-pts(1,1)+1, pts(3,2)-pts(3,1)+1, 1), &
                         coorZ(pts(1,2)-pts(1,1)+1, pts(3,2)-pts(3,1)+1, 1))
                    nSymm = dims(1)*dims(3)
                 else if (pts(3,1) == pts(3, 2)) then 
                    allocate(&
                         coorX(pts(1,2)-pts(1,1)+1, pts(2,2)-pts(2,1)+1, 1), &
                         coorY(pts(1,2)-pts(1,1)+1, pts(2,2)-pts(2,1)+1, 1), &
                         coorZ(pts(1,2)-pts(1,1)+1, pts(2,2)-pts(2,1)+1, 1))
                    nSymm = dims(1)*dims(2)
                 end if

                 call cg_coord_read_f(cg, base, iZone, "CoordinateX", RealDouble, pts(:,1), pts(:,2), coorX, ierr)
                 if (ierr .eq. CG_ERROR) call cg_error_exit_f
                 call cg_coord_read_f(cg, base, iZone, "CoordinateY", RealDouble, pts(:,1), pts(:,2), coorY, ierr)
                 if (ierr .eq. CG_ERROR) call cg_error_exit_f
                 call cg_coord_read_f(cg, base, iZone, "CoordinateZ", RealDouble, pts(:,1), pts(:,2), coorZ, ierr)
                 if (ierr .eq. CG_ERROR) call cg_error_exit_f

                 ! Determine what the axis of this symmetry plane
                 ! is. We look at the two end points of the range
                 symmSum = zero
                 symmSum(1) = sum(abs(coorX))
                 symmSum(2) = sum(abs(coorY))
                 symmSum(3) = sum(abs(coorZ))

                 if (symmSum(1) / nSymm < symmTol) then 
                    tmpSym = 1
                 else if (symmSum(2) / nSymm < symmTol) then 
                    tmpSym = 2
                 else if (symmSum(3) / nSymm < symmTol) then
                    tmpSym = 3
                 end if
                 
                 if (iSymm == 0) then ! Currently un-assigned:
                    iSymm = tmpSym
                 end if

                 if (tmpSym /= iSymm) then
                    print *, 'Error: detected more than 1 symmetry plane direction.'
                    print *, 'pyWarpUnstruct cannot handle this'
                    stop
                 end if
                 deallocate(coorX, coorY, coorZ)
              end if
           end if
        end do bocoLoop_one
     end do zoneLoop_one

     ! Allocate space for storage of the surface nodes and the
     ! connectivity
     allocate(wallPoints(3*nWall))
     allocate(wallConn(4*nConn))

     ! Reset the counters here
     nWall = 0
     offSet = 0
     nConn = 0
     nodeCount = 0
     ii = 0
     zoneLoop_two: do iZone=1,nZones
        dims(1) = zones(iZone)%il
        dims(2) = zones(iZone)%jl
        dims(3) = zones(iZone)%kl

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

        bocoLoop_two: do boco=1, size(zones(iZone)%bocos)
           bocoType = zones(iZone)%bocos(boco)%type
           pts = zones(iZone)%bocos(boco)%ptRange
           if (BCTypeName(bocoType) == 'BCWallViscous' .or. &
                BCTypeName(bocoType) == 'BCWallInviscid' .or. &
                BCTypeName(bocoType) == 'BCWall' .or. &
                BCTypeName(bocoType) == 'BCWallViscousHeatFlux' .or. &
                BCTypeName(bocoType) == 'BCWallViscousIsothermal') then

              ! Flag the wall nodes for walls:
              do k=pts(3, 1), pts(3, 2)
                 do j=pts(2, 1), pts(2, 2)
                    do i=pts(1, 1), pts(1, 2)
                       ! Wall nodes in global ordering
                       wallNodes(offset + (k-1)*dims(1)*dims(2) + &
                                          (j-1)*dims(1) + &
                                          i) = 1

                       ! List of all wall nodes
                       wallPoints(3*nWall+1) = coorX(i, j, k)
                       wallPoints(3*nWall+2) = coorY(i, j, k)
                       wallPoints(3*nWall+3) = coorZ(i, j, k)
                       nWall = nWall + 1
                    end do
                 end do
              end do

              lowFace = .False.
              ! Determine generic face size
              if (pts(1,1) == pts(1,2)) then ! iMin/iMax
                 ni = pts(2,2) - pts(2,1) + 1
                 nj = pts(3,2) - pts(3,1) + 1
                 if (pts(1,1) == 1) then 
                    lowFace = .True.
                 end if

              else if (pts(2,1) ==pts(2,2)) then !jMin/jMax
                 ni = pts(1,2) - pts(1,1) + 1
                 nj = pts(3,2) - pts(3,1) + 1
                 if (pts(2,1) == 1) then 
                    lowFace = .True.
                 end if

              else ! kMin/kMax
                 ni = pts(1,2) - pts(1,1) + 1
                 nj = pts(2,2) - pts(2,1) + 1

                 if (pts(3,1) == 1) then 
                    lowFace = .True.
                 end if

              end if
          
              ! Loop over generic face size...We are doing 1 based
              ! ordering. If it is low face normal ordering: 
              !
              ! i, j+1 +-----+ i+1, j+1
              !   n4   |     | n3
              !        +-----+
              !       i,j    i+1, j
              !       n1     n2
              !
              ! Otherwise, we flip the ordering to be n1, n4, n3, n2
              if (lowFace) then 
                 do j=0,nj-2
                    do i=0,ni-2
                       wallConn(4*nConn+1) = nodeCount + (j  )*ni + i + 1     ! n1
                       wallConn(4*nConn+2) = nodeCount + (j  )*ni + i + 1 + 1 ! n2
                       wallConn(4*nConn+3) = nodeCount + (j+1)*ni + i + 1 + 1 ! n3
                       wallConn(4*nConn+4) = nodeCount + (j+1)*ni + i + 1     ! n4
                       nConn = nConn + 1
                    end do
                 end do
              else ! Flip the orientation:
                 do j=0,nj-2
                    do i=0,ni-2
                       wallConn(4*nConn+1) = nodeCount + (j  )*ni + i + 1     ! n1
                       wallConn(4*nConn+2) = nodeCount + (j+1)*ni + i + 1     ! n4
                       wallConn(4*nConn+3) = nodeCount + (j+1)*ni + i + 1 + 1 ! n3
                       wallConn(4*nConn+4) = nodeCount + (j  )*ni + i + 1 + 1 ! n2
                       nConn = nConn + 1
                    end do
                 end do
              end if
              nodeCount = nodeCount + ni*nj
           end if

           if (BCTypeName(bocoType) == 'BCSymmetryPlane') then 
              ! Hard zero the sym plane in case our mesh is bad. We
              ! already know what iSymm is from above.
              do k=pts(3, 1), pts(3, 2)
                 do j=pts(2, 1), pts(2, 2)
                    do i=pts(1, 1), pts(1, 2)
                       allNodes(isymm, offset + (k-1)*dims(1)*dims(2) + (j-1)*dims(1) + i) = zero
                    end do
                 end do
              end do
           end if
        end do bocoLoop_two

        offset = offset + dims(1)*dims(2)*dims(3)

        ! And free the temporary arrays
        deallocate(coorX, coorY, coorZ)

     end do zoneLoop_two

     call cg_close_f(cg, ierr)
     if (ierr .eq. CG_ERROR) call cg_error_exit_f
     deallocate(sizes)
  else
     ! Allocate these to zero so we can just blindly dealloc later
     allocate(wallPoints(0), wallConn(0))
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
