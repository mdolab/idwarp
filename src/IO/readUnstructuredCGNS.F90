subroutine readUnstructuredCGNS(cgns_file)

  use precision
  use communication
  use gridData
  use gridInput  
  use unStructuredCGNSGrid

  implicit none
  include 'cgnslib_f.h'

  ! Input Arguments
  character*(*),intent(in) :: cgns_file

  ! CGNS Variabls
  integer(kind=intType) :: i, ii, istart, iend, localsize, iProc, iZone, offset
  integer(kind=intType):: ierr, base, dims(3), cg!, cursize, coorsize
  integer(kind=intType):: nNodes, nCells
  integer(kind=intType):: CellDim, PhysDim
  integer(kind=intType) :: nbases, start(3)!, tmpSym, nSymm
  character*32 :: zoneName, bocoName, famName, connectName, donorName, basename
  character*32 :: secName

  integer(kind=intType) :: nbocos, boco, index
  integer(kind=intType) :: nVertices,nElements
  integer(kind=intType) :: zoneType,dataType,sec,type
  integer(kind=intType) :: nCoords,nSections,nElem,nConn
  integer(kind=intType) :: eBeg,eEnd,nBdry, parentFlag
  integer(kind=intType) :: bocoType,ptsettype,nbcelem
  integer(kind=intType) :: normalIndex,normalListFlag
  integer(kind=intType) :: normDataType,nDataSet
  real(kind=realType), dimension(:), allocatable :: coorX, coorY, coorZ
#ifdef USE_COMPLEX
  complex(kind=realType), dimension(:, :), allocatable :: allNodes, localNodes
#else
  real(kind=realType), dimension(:, :), allocatable :: allNodes, localNodes
#endif
  integer(kind=intType), dimension(:), allocatable :: wallNodes, localWallNodes
  integer(kind=intType), dimension(:, :), allocatable :: sizes

  integer(kind=intType) :: status(MPI_STATUS_SIZE)
  integer(kind=intType):: surfSecCounter,volSecCounter, famID

  iSymm = 0
  ! ---------------------------------------
  !           Open CGNS File
  ! ---------------------------------------
  if (myID == 0) then
     print *, ' -> Reading CGNS File: ', cgns_file

     call cg_open_f(cgns_file, CG_MODE_READ, cg, ierr)
     if (ierr .eq. CG_ERROR) call cg_error_exit_f

     call cg_nbases_f(cg, nbases, ierr)
     if (ierr .eq. CG_ERROR) call cg_error_exit_f

     if (nbases .gt. 1) then
        print *, ' ** Warning: pyWarpUstruct only reads the first base in a cgns file'
     end if

     base = 1_intType

     !       *** base attribute:  GOTO base node
     call cg_goto_f(cg, base, ierr, 'end')
     if (ierr .eq. CG_ERROR) call cg_error_exit_f

     ! Check the cell and physical dimensions of the bases. 
     ! Both should be 3. 
     call cg_base_read_f(cg,base,baseName,cellDim,physDim,ierr)
     if (ierr .eq. ERROR) call cg_error_exit_f

     if( cellDim .ne.3 .or. physDim .ne. 3) then
        print *,' ** Warning:  pyWarpUstruct only accepts 3-d data'
        stop
     end if

     ! read the number of zones in the file
     call cg_nzones_f(cg, base, nzones, ierr)
     if (ierr .eq. ERROR) call cg_error_exit_f

     print *, '   -> Number of Zones:', nzones

     ! allocate the zonal storage
     allocate(zones(nZones))

     ! do a check that we have all unstructured zones
     do iZone = 1,nZones
        ! Check the zone type. This should be Unstructured. 
        call cg_zone_type_f(cg,base, iZone, zoneType,ierr)
        if (ierr .eq. ERROR) call cg_error_exit_f

        if (zoneType .ne. Unstructured) then
           print *,"Unstructured zone expected...stopping"
           stop
        end if
     end do

     ! Now count up the total number of nodes
     nNodes = 0
     nCells = 0
     allocate(sizes(3, nZones))

     do iZone=1, nZones
        call cg_zone_read_f(cg, base, iZone, zoneName, dims, ierr)
        nVertices = dims(1)
        nElements = dims(2)
        zones(iZone)%nVertices = nVertices
        zones(iZone)%nElements = nElements
        zones(iZone)%zoneName = zoneName
        sizes(:, iZone) = dims(1:3)
        nNodes = nNodes + dims(1)
        nCells = nCells + dims(2)
     end do

     ! Now we know the total number of nodes we can allocate the final
     ! required space and read them in. 
     allocate(allNodes(3, nNodes))
     allocate(wallNodes(nNodes))
     wallNodes = 0

     ii = 0
     familyList(:) = ''
     nWallFamilies = 0
     offset = 0

     ! loop over the zones and read the nodes
     zoneLoop: do iZone = 1,nZones
        call cg_zone_read_f(cg, base, iZone, zoneName, dims, ierr)
        if (ierr .eq. CG_ERROR) call cg_error_exit_f
        allocate(coorX(dims(1)), coorY(dims(1)), coorZ(dims(1)))

        ! Read the x,y,z-coordinates.
        ! This assumes the coords are in double precision.
        ! Note that CGNS starts the numbering at 
        ! 1 even if C is used. 

        call cg_coord_read_f(cg,base,iZone,"CoordinateX",&
             & realDouble, 1,dims(1), coorX,ierr)
        if (ierr .eq. CG_ERROR) call cg_error_exit_f
        
        call cg_coord_read_f(cg,base,iZone,"CoordinateY",&
             & realDouble, 1,dims(1), coorY,ierr)
        if (ierr .eq. CG_ERROR) call cg_error_exit_f
           
        call cg_coord_read_f(cg,base,iZone,"CoordinateZ",&
             & realDouble, 1,dims(1), coorZ,ierr)
        if (ierr .eq. CG_ERROR) call cg_error_exit_f

        ! Now stack all of the zones in one array
        do i=1, dims(1)
           ii = ii + 1
#ifdef USE_COMPLEX
           allNodes(1, ii) = cmplx(coorX(i), 0.0)
           allNodes(2, ii) = cmplx(coorY(i), 0.0)
           allNodes(3, ii) = cmplx(coorZ(i), 0.0)
#else
           allNodes(1, ii) = coorX(i)
           allNodes(2, ii) = coorY(i)
           allNodes(3, ii) = coorZ(i)
#endif
        end do

        ! Now loop over the sections in the unstructured file and
        ! figure out the boundary nodes.

        ! Determine the number of sections for this zone. Note that 
        ! surface elements can be stored in a volume zone, but they 
        ! are NOT taken into account in the number obtained from 
        ! cg_zone_read. 

        call cg_nsections_f(cg, base, iZone, nSections, ierr)
        if (ierr .eq. CG_ERROR) call cg_error_exit_f
        print *, '   -> nSections',nSections

        zones(iZone)%nSections = nSections

        ! allocate the logical to track the section type
        allocate(zones(iZone)%isVolumeSection(nSections))
        zones(iZone)%isVolumeSection = .True.

        ! loop over the number of sections and determine which
        ! sections are volume elements and which are surface elements

        zones(iZone)%nVolSec = 0
        zones(iZone)%nSurfSec = 0

        do sec=1,nSections
           ! Determine the element type and set the pointer for the 
           ! connectivity accordingly. 
           call cg_section_read_f(cg, base, iZone, sec, secName, &
                type, eBeg, eEnd, nBdry, parentFlag, ierr)
           if (ierr .eq. CG_ERROR) call cg_error_exit_f

           select case (type)
           case (TETRA_4)
              zones(iZone)%nVolSec = zones(iZone)%nVolSec+1
              zones(iZone)%isVolumeSection(sec) = .true.
           case (PYRA_5)
              zones(iZone)%nVolSec = zones(iZone)%nVolSec+1
              zones(iZone)%isVolumeSection(sec) = .true.
           case (PENTA_6)
              zones(iZone)%nVolSec = zones(iZone)%nVolSec+1
              zones(iZone)%isVolumeSection(sec) = .true.
           case (HEXA_8)
              zones(iZone)%nVolSec = zones(iZone)%nVolSec+1
              zones(iZone)%isVolumeSection(sec) = .true.
           case (TRI_3)
              zones(iZone)%nSurfSec = zones(iZone)%nSurfSec+1
              zones(iZone)%isVolumeSection(sec) = .false.
           case (QUAD_4)
              zones(iZone)%nSurfSec =  zones(iZone)%nSurfSec+1
              zones(iZone)%isVolumeSection(sec) = .false.
           case  default
              if (myID == 0) then
                 print *, "Unsupported element encountered....exiting", type
                 stop
              end if
           end select
        end do

        allocate(zones(iZone)%volumeSections(zones(iZone)%nVolSec))
        allocate(zones(iZone)%surfaceSections(zones(iZone)%nSurfSec))

        ! Loop back over the number of sections and read the element 
        ! connectivities. As CGNS starts the numbering at 1 the 
        ! for-loop starts at 1 as well. 

        volSecCounter = 1
        surfSecCounter = 1

        do sec=1, nSections
           ! Determine the element type and set the pointer for the 
           ! connectivity accordingly. 
           call cg_section_read_f(cg, base, iZone, sec, secName, &
                type, eBeg, eEnd, nBdry, parentFlag, ierr)
           if (ierr .eq. CG_ERROR) call cg_error_exit_f
           nElem = eEnd - eBeg + 1

           ! Get the nodes per element "npe"
           call cg_npe_f(type, nconn, ierr)
           if (ierr .eq. CG_ERROR) call cg_error_exit_f

           if (zones(iZone)%isVolumeSection(sec)) then
              ! this is a group of volume elements
              zones(iZone)%volumeSections(volSecCounter)%nElem = nElem
              zones(iZone)%volumeSections(volSecCounter)%secName = secName
              zones(iZone)%volumeSections(volSecCounter)%elemType = type
              zones(iZone)%volumeSections(volSecCounter)%nConn = nConn 
              zones(iZone)%volumeSections(volSecCounter)%elemStart = eBeg
              zones(iZone)%volumeSections(volSecCounter)%elemEnd = eEnd
              volSecCounter = volSecCounter+1
           else
              ! this is a group of surface elements
              zones(iZone)%surfaceSections(surfSecCounter)%nElem = nElem
              zones(iZone)%surfaceSections(surfSecCounter)%secName = secName
              zones(iZone)%surfaceSections(surfSecCounter)%elemType = type
              zones(iZone)%surfaceSections(surfSecCounter)%nConn = nConn 
              zones(iZone)%surfaceSections(surfSecCounter)%elemStart = eBeg
              zones(iZone)%surfaceSections(surfSecCounter)%elemEnd = eEnd   
              surfSecCounter = surfSecCounter + 1
           end if
        end do

        !determine the number of boundary conditions for this zone. 
        call cg_nbocos_f(cg, base, iZone, nBocos, ierr)
        if (ierr .eq. CG_ERROR) call cg_error_exit_f

        ! Loop over the number of boundary conditions.
        zones(iZone)%nBocos = nBocos
        allocate(zones(iZone)%BCInfo(nBocos),STAT=ierr)

        bocoLoop: do boco=1,nBocos
           ! Read the info for this boundary condition. 
           call cg_boco_info_f(cg, base, iZone, boco, boconame, bocotype, &
                ptsettype, nbcElem, NormalIndex, NormalListFlag, datatype, ndataset,&
                ierr)
           if (ierr .eq. CG_ERROR) call cg_error_exit_f

           zones(iZone)%BCInfo(boco)%nBCElem = nBCElem
           zones(iZone)%BCInfo(boco)%BCName =  bocoName
           zones(iZone)%BCInfo(boco)%BCType = bocoType
           allocate(zones(iZone)%BCInfo(boco)%BCElements(nBCElem))

           ! Read the element ID's.
           call cg_boco_read_f(cg, base, iZone, boco,&
                zones(iZone)%BCInfo(boco)%BCElements, NULL, ierr)

           if (ierr .eq. CG_ERROR) call cg_error_exit_f

           ! Now that we have the boundary condition info, read the family name
           ! to associate it with a specific section

           call cg_goto_f(cg, base, ierr, "Zone_t", iZone, "ZoneBC_t", 1, "BC_t",&
                boco, "end")
           if (ierr .eq. CG_ERROR) call cg_error_exit_f

           if (ierr == 0) then ! Node exists
              call cg_famname_read_f(famName, ierr)
              if (ierr .eq. CG_ERROR) call cg_error_exit_f
              if (ierr == 0) then
                 !search for matching section
                 ! Loop over the number of sections and compare the family
                 ! name to the section name. If they match, we have found
                 ! the section that this BC belongs to.

                 volSecCounter = 1
                 surfSecCounter = 1
                 do sec=1,nSections
                    if (.not. zones(iZone)%isVolumeSection(sec)) then
                       ! we have a surface section
                       zones(iZone)%BCInfo(boco)%famName = famName
 
                      if(zones(iZone)%surfaceSections(surfSecCounter)%secName==&
                            famName)then
                          ! this is the correct section
                          ! Set BCType and Family Name
                          zones(iZone)%surfaceSections(surfSecCounter)%BCFamily = famName
                          zones(iZone)%surfaceSections(surfSecCounter)%BCType = bocoType
                       end if
                       surfSecCounter = surfSecCounter + 1
                    end if
                 end do
              end if
           end if
        end do bocoLoop

        ! And free the temporary arrays
        deallocate(coorX, coorY, coorZ)
     end do zoneLoop

     ! We can also generate the wall family list
     nwallFamilies = 0
     familyList(:) = ''

     do iZone=1, nzones
        nSections = zones(iZone)%nSections
        surfSecCounter = 1
        do sec=1, nSections
           if (.not. zones(iZone)%isVolumeSection(sec)) then
              bocoType = zones(iZone)%surfaceSections(surfSecCounter)%BCType
              FamName = zones(iZone)%surfaceSections(surfSecCounter)%BCFamily

              zones(iZone)%surfaceSections(surfSecCounter)%isWallBC=.False.
              zones(iZone)%surfaceSections(surfSecCounter)%isBoundaryBC=.False.
              zones(iZone)%surfaceSections(surfSecCounter)%isSymmBC=.False.

              if (BCTypeName(bocoType) == 'BCWallViscous' .or. &
                   BCTypeName(bocoType) == 'BCWallInviscid' .or. &
                   BCTypeName(bocoType) == 'BCWall' .or. &
                   BCTypeName(bocoType) == 'BCWallViscousHeatFlux' .or. &
                   BCTypeName(bocoType) == 'BCWallViscousIsothermal') then
                 zones(iZone)%surfaceSections(surfSecCounter)%isWallBC=.True.

                 !this is a wall family, add to the family list
                 call checkInFamilyList(familyList, famName, nwallFamilies, famID)
                 if (famID == 0) then
                    nwallFamilies = nwallFamilies + 1
                    familyList(nwallFamilies) = famName
                 end if
              elseif(BCTypeName(bocotype) == 'BCFarfield')  then
                 zones(iZone)%surfaceSections(surfSecCounter)%isBoundaryBC=.True.
              elseif(BCTypeName(bocotype) == 'BCSymmetryPlane') then
                 zones(iZone)%surfaceSections(surfSecCounter)%isSymmBC = .True.
              else
                 print *,'Unrecongnized boundary Type... exiting: ',BCTypeName(bocotype)
                 stop
              end if
              surfSecCounter = surfSecCounter + 1
           end if
        end do
     end do
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

     ! calculate the start and end indices for the root process
     istart = 0 + 1
     iend   = int(nNodes*(one/nProc))
     localSize = iend - istart + 1

     ! Allocate the local storage for the nodes
     allocate(localNodes(3, localSize))
     allocate(localWallNodes(localSize))

     ! Just copy for the root proc:
     localNodes(:, :) = allNodes(:, 1:localSize)
     localWallNodes(:) = wallNodes(1:localSize)

     ! Loop over the remainer of the procs and send the nodes they need
     do iProc=1, nProc-1
        ! compute the node range for the target processor
        istart = int(nNodes*(dble(iProc)/nProc)) +1 
        iend   = int(nNodes*(dble(iProc+1)/nProc))
        localSize = iend - istart + 1

        ! Post the send message for the nodes
#ifdef USE_COMPLEX
        call MPI_Send(allNodes(:, iStart), localSize*3, MPI_COMPLEX16, iProc, &
             11, warp_comm_world, ierr)
#else
        call MPI_Send(allNodes(:, iStart), localSize*3, MPI_REAL8, iProc, &
             11, warp_comm_world, ierr)
#endif
        call EChk(ierr, __FILE__, __LINE__)

        ! post the send message for the wall node map.
        call MPI_Send(wallNodes(iStart), localSize, MPI_INTEGER4, iProc, &
             12, warp_comm_world, ierr)
        call EChk(ierr, __FILE__, __LINE__)
     end do
     deallocate(allNodes)
     deallocate(wallNodes)
  else

     ! compute the incoming node range for this processor
     istart = int(nNodes*(dble(myid)/nProc)) +1 
     iend   = int(nNodes*(dble(myid+1)/nProc))
     localSize = iend - istart + 1

     ! allocate the storage for the incoming data
     allocate(localNodes(3, localSize))
     allocate(localWallNodes(localSize))

     ! Post the receive on each processor
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

end subroutine readUnstructuredCGNS


