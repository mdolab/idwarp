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


  !   integer(kind=intType):: ierr
  !   ! integer(kind=intType):: index_base, index_zone
  !   ! integer(kind=intType), dimension(3,1):: isize
  !   ! integer(kind=intType), dimension(3):: irmin, irmax

  !   character*32 :: zoneName,name,secName,bocoName,famName
  !   !character*32 :: baseName
  ! CGNS Variabls
  integer(kind=intType) :: i, ii, istart, iend, localsize, iProc, iZone, offset
  integer(kind=intType):: ierr, base, dims(3), cg!, cursize, coorsize
  integer(kind=intType):: nNodes, nCells
  integer(kind=intType):: CellDim, PhysDim
  integer(kind=intType) :: nbases, start(3)!, tmpSym, nSymm
  character*32 :: zoneName, bocoName, famName, connectName, donorName, basename
  character*32 :: secName
  !   integer(kind=intType) :: fileIndex, base, nbases
  integer(kind=intType) :: nbocos, boco, index
  integer(kind=intType) :: nVertices,nElements
  integer(kind=intType) :: zoneType,dataType,sec,type
  integer(kind=intType) :: nCoords,nSections,nElem,nConn
  integer(kind=intType) :: eBeg,eEnd,nBdry, parentFlag
  integer(kind=intType) :: bocoType,ptsettype,nbcelem
  integer(kind=intType) :: normalIndex,normalListFlag
  integer(kind=intType) :: normDataType,nDataSet
  !integer(kind=intType) :: i,ii!,j,counter
  !   !integer(kind=8), dimension(3)::sizes
  !   !integer(kind=intType), dimension(3)::sizes
  !   integer(kind=intType), dimension(:,:), allocatable::conn
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
  !   !integer(kind=intType), dimension(:), allocatable:: BCElements
  !   !integer, dimension(:), allocatable::conn,parentData

  iSymm = 0
  ! ---------------------------------------
  !           Open CGNS File
  ! ---------------------------------------
  if (myID == 0) then
     print *, ' -> Reading CGNS File: ', cgns_file


     ! if (CG_BUILD_64BIT .ne. 0) then
     !    print *,'will not work in 64-bit mode',CG_BUILD_64BIT 
     !    stop
     ! endif

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

        if(zoneType .ne. Unstructured)then
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
        !print *,'Zone Sizes',zoneName, dims
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

     ii=0
     familyList(:) = ''
     nWallFamilies = 0
     offset = 0
     ! loop over the zones and read the nodes
     zoneLoop: do iZone = 1,nZones
        call cg_zone_read_f(cg, base, iZone, zoneName, dims, ierr)
        if (ierr .eq. CG_ERROR) call cg_error_exit_f
        ! allocate(gridDoms(zone)%points(nVertices,physDim),STAT=ierr)
        ! allocate(gridDoms(zone)%points0(nVertices,physDim),STAT=ierr)
        ! allocate(gridDoms(zone)%dx(nVertices,physDim),STAT=ierr)
        ! allocate(gridDoms(zone)%isSurfaceNode(nVertices),STAT=ierr)
        allocate(coorX(dims(1)), &
             coorY(dims(1)), &
             coorZ(dims(1)))

        ! Read the x,y,z-coordinates.
        ! This assumes the coords are in double precision.
        ! Note that CGNS starts the numbering at 
        ! 1 even if C is used. 

        call cg_coord_read_f(cg,base,iZone,"CoordinateX",&
             & realDouble, 1,dims(1), coorX,ierr)

        call cg_coord_read_f(cg,base,iZone,"CoordinateY",&
             & realDouble, 1,dims(1), coorY,ierr)

        call cg_coord_read_f(cg,base,iZone,"CoordinateZ",&
             & realDouble, 1,dims(1), coorZ,ierr)

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

        call cg_nsections_f(cg,base,iZone,nSections,ierr)
        print *, '   -> nSections',nSections

        zones(iZone)%nSections = nSections
        ! allocate the logical to track the section type
        allocate(zones(iZone)%isVolumeSection(nSections),STAT=ierr)
        zones(iZone)%isVolumeSection=.True.

        ! loop over the number of sections and determine which
        ! sections are volume elements and which are surface elements

        zones(iZone)%nVolSec = 0
        zones(iZone)%nSurfSec = 0

        do sec=1,nSections
           ! Determine the element type and set the pointer for the 
           ! connectivity accordingly. 
           call cg_section_read_f(cg,base,iZone,sec,secName,&
                type,eBeg,eEnd,nBdry,parentFlag,ierr)

           select case (type)

           case (TETRA_4)
              zones(iZone)%nVolSec = zones(iZone)%nVolSec+1
              zones(iZone)%isVolumeSection(sec)=.true.
           case (PYRA_5)
              zones(iZone)%nVolSec = zones(iZone)%nVolSec+1
              zones(iZone)%isVolumeSection(sec)=.true.
           case (PENTA_6)
              zones(iZone)%nVolSec = zones(iZone)%nVolSec+1
              zones(iZone)%isVolumeSection(sec)=.true.
           case (HEXA_8)
              zones(iZone)%nVolSec = zones(iZone)%nVolSec+1
              zones(iZone)%isVolumeSection(sec)=.true.
           case (TRI_3)
              zones(iZone)%nSurfSec = zones(iZone)%nSurfSec+1
              zones(iZone)%isVolumeSection(sec)=.false.
           case (QUAD_4)
              zones(iZone)%nSurfSec =  zones(iZone)%nSurfSec+1
              zones(iZone)%isVolumeSection(sec)=.false.
           case  default
              if (myID == 0) then
                 print *, "Unsupported element encountered....exiting", type
                 stop
              end if
           end select

        end do

        allocate(zones(iZone)%volumeSections(zones(iZone)%nVolSec),STAT=ierr)
        allocate(zones(iZone)%surfaceSections(zones(iZone)%nSurfSec),STAT=ierr)

        ! Loop over the number of sections and read the element 
        ! connectivities. As CGNS starts the numbering at 1 the 
        ! for-loop starts at 1 as well. 

        volSecCounter = 1
        surfSecCounter = 1

        do sec=1,nSections
           ! Determine the element type and set the pointer for the 
           ! connectivity accordingly. 
           call cg_section_read_f(cg,base,iZone,sec,secName,&
                type,eBeg,eEnd,nBdry,parentFlag,ierr)
           nElem = eEnd-eBeg+1
           !print *,'type',sec,secName,type,eBeg,eEnd,nBdry,parentFlag
           !print *,'Types',HEXA_8,QUAD_4

           call cg_npe_f(type, nconn, ierr)

           if (zones(iZone)%isVolumeSection(sec)) then
              ! this is a group of volume elements
              ! allocate(zones(iZone)%volumeSections(volSecCounter)%elemPtr(nElem+1),STAT=ierr)
              ! allocate(zones(iZone)%volumeSections(volSecCounter)%elemConn(nConn*nElem),STAT=ierr)
              zones(iZone)%volumeSections(volSecCounter)%nElem = nElem
              zones(iZone)%volumeSections(volSecCounter)%secName = secName
              zones(iZone)%volumeSections(volSecCounter)%elemType = type
              zones(iZone)%volumeSections(volSecCounter)%nConn = nConn 
              zones(iZone)%volumeSections(volSecCounter)%elemStart = eBeg
              zones(iZone)%volumeSections(volSecCounter)%elemEnd = eEnd
           else
              ! this is a group of surface elements
              ! allocate(zones(iZone)%surfaceSections(surfSecCounter)%elemPtr(nElem+1),STAT=ierr)
              ! allocate(zones(iZone)%surfaceSections(surfSecCounter)%elemConn(nConn*nElem),STAT=ierr)

              ! ! allocate(zones(iZone)%surfaceSections(surfSecCounter)%elements(nConn,nElem),STAT=ierr)
              ! allocate(zones(iZone)%surfaceSections(surfSecCounter)%elemCenter(nElem,physDim),STAT=ierr)
              ! allocate(zones(iZone)%surfaceSections(surfSecCounter)%elemArea(nElem,physDim),STAT=ierr)
              ! allocate(zones(iZone)%surfaceSections(surfSecCounter)%elemNormal(nElem,physDim),STAT=ierr)
              ! allocate(zones(iZone)%surfaceSections(surfSecCounter)%elemAreaMag(nElem),STAT=ierr)
              zones(iZone)%surfaceSections(surfSecCounter)%nElem = nElem
              zones(iZone)%surfaceSections(surfSecCounter)%secName = secName
              zones(iZone)%surfaceSections(surfSecCounter)%elemType = type
              zones(iZone)%surfaceSections(surfSecCounter)%nConn = nConn 
              zones(iZone)%surfaceSections(surfSecCounter)%elemStart = eBeg
              zones(iZone)%surfaceSections(surfSecCounter)%elemEnd = eEnd   

           end if

           !         ! call cg_npe_f(type, nconn, ierr)
           !         allocate(conn(nConn,nElem),STAT=ierr)
           !         conn = 0
           !         ! allocate(parentData(nElem,nConn),STAT=ierr)
           !         ! Read the connectivity. Again, the node numbering of the 
           !         ! connectivities start at 1. If internally a starting
           !         ! index of 0 is used (typical for C-codes) 1 must be
           !         ! substracted from the connectivities read. 

           !         call cg_elements_read_f(fileIndex,base,iZone,sec,conn,parentFlag,ierr)
           !         ! if (myID == 0) then
           !         !    print *, 'conn',shape(conn),conn(:,1)
           !         ! end if

           if (zones(iZone)%isVolumeSection(sec)) then
              !            counter = 1
              !            do i = 1,nElem
              !               zones(iZone)%volumeSections(volSecCounter)%elemPtr(i)=counter
              !               do j = 1,nConn
              !                  zones(iZone)%volumeSections(volSecCounter)%elemConn(counter)=conn(j,i)
              !                  counter = counter +1
              !                  !zones(iZone)%volumeSections(volSecCounter)%elements = conn
              !               end do
              !            end do
              !            zones(iZone)%volumeSections(volSecCounter)%elemPtr(i)=counter
              volSecCounter = volSecCounter+1
           else
              !            !zones(iZone)%surfaceSections(surfSecCounter)%elements = conn
              !            counter = 1
              !            do i = 1,nElem
              !               zones(iZone)%surfaceSections(surfSecCounter)%elemPtr(i)=counter
              !               do j = 1,nConn
              !                  zones(iZone)%surfaceSections(surfSecCounter)%elemConn(counter)=conn(j,i)
              !                  counter = counter +1
              !               end do
              !               !print *,'counter',i,counter,zones(iZone)%surfaceSections(surfSecCounter)%elemPtr(i)
              !            end do
              !            !print *,'end',i,counter
              !            zones(iZone)%surfaceSections(surfSecCounter)%elemPtr(i)=counter
              surfSecCounter = surfSecCounter + 1
           end if
           !         deallocate(conn,STAT=ierr)
        end do

        !determine the number of boundary conditions for this zone. 
        call cg_nbocos_f(cg,base,iZone,nBocos,ierr)

        ! Loop over the number of boundary conditions.
        zones(iZone)%nBocos = nBocos
        allocate(zones(iZone)%BCInfo(nBocos),STAT=ierr)

        bocoLoop: do boco=1,nBocos
           ! Read the info for this boundary condition. 
           call cg_boco_info_f(cg, base, iZone, boco, boconame, bocotype, &
                ptsettype, nbcElem, NormalIndex, NormalListFlag, datatype, ndataset,&
                ierr)
           if (ierr .eq. CG_ERROR) call cg_error_exit_f

           !print *,'boco',boco,bocoType,BCTypeName(bocoType)
           zones(iZone)%BCInfo(boco)%nBCElem = nBCElem
           zones(iZone)%BCInfo(boco)%BCName =  bocoName
           zones(iZone)%BCInfo(boco)%BCType = bocoType
           allocate(zones(iZone)%BCInfo(boco)%BCElements(nBCElem),STAT=ierr)

           ! Read the element ID's.
           call cg_boco_read_f(cg,base,iZone,boco,&
                zones(iZone)%BCInfo(boco)%BCElements,NULL,ierr)

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
                       !print *,'names',famName,zones(iZone)%surfaceSections(surfSecCounter)%secName
                       if(zones(iZone)%surfaceSections(surfSecCounter)%secName==&
                            famName)then
                          !print *,'names Match',famName,zones(iZone)%surfaceSections(surfSecCounter)%secName,surfSecCounter,bocoType
                          ! this is the correct section
                          ! Set BCType and Family Name
                          zones(iZone)%surfaceSections(surfSecCounter)%BCFamily = famName
                          zones(iZone)%surfaceSections(surfSecCounter)%BCType = bocoType
                          ! else
                          !    ! this is not the correct section
                       end if
                       surfSecCounter = surfSecCounter + 1
                    end if
                 end do

                 !print *,'BCs',iZone,index,zones(iZone)%surfaceSections(index)%BCFamily,zones(iZone)%surfaceSections(index)%BCType
                 ! else
                 !    !family not present set logical
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
        do sec=1,nSections
           if (.not. zones(iZone)%isVolumeSection(sec)) then
              bocoType = zones(iZone)%surfaceSections(surfSecCounter)%BCType
              FamName = zones(iZone)%surfaceSections(surfSecCounter)%BCFamily
              zones(iZone)%surfaceSections(surfSecCounter)%isWallBC=.False.
              zones(iZone)%surfaceSections(surfSecCounter)%isBoundaryBC=.False.
              zones(iZone)%surfaceSections(surfSecCounter)%isSymmBC=.False.
              !print *,'bocotype',bocoType,BCTypeName(bocoType),surfSecCOunter
              if (BCTypeName(bocoType) == 'BCWallViscous' .or. &
                   BCTypeName(bocoType) == 'BCWallInviscid' .or. &
                   BCTypeName(bocoType) == 'BCWall' .or. &
                   BCTypeName(bocoType) == 'BCWallViscousHeatFlux' .or. &
                   BCTypeName(bocoType) == 'BCWallViscousIsothermal') then
                 !print *,'is Wall',iZone,surfSecCounter
                 zones(iZone)%surfaceSections(surfSecCounter)%isWallBC=.True.
                 !this is a wall family, add to the family list
                 call checkInFamilyList(familyList, famName, nwallFamilies, famID)
                 if (famID == 0) then
                    nwallFamilies = nwallFamilies + 1
                    familyList(nwallFamilies) = famName
                 end if
              elseif(BCTypeName(bocotype) == 'BCFarfield')  then
                 !print *,'is farfield',iZone,surfSecCounter
                 zones(iZone)%surfaceSections(surfSecCounter)%isBoundaryBC=.True.
              elseif(BCTypeName(bocotype) == 'BCSymmetryPlane') then
                 !print *,' is symm',iZone,surfSecCounter
                 zones(iZone)%surfaceSections(surfSecCounter)%isSymmBC = .True.
                 hasSymmetry = .True.
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


  !   ! initialize the Symmetry condition to false
  !   hasSymmetry = .False.

  !   ! Allocate the GridDoms Data structure
  !  ! allocate(zones(nZones),STAT=ierr)
  !   call allocateGridDoms(nZones)

  !   ! loop over the zones and read
  !   do zone = 1,nZones

  !      ! Check the zone type. This should be Unstructured. 
  !      call cg_zone_type_f(fileIndex,base, zone, zoneType,ierr)
  !      if (ierr .eq. ERROR) call cg_error_exit_f
  !      if (myID == 0) then
  !         print *, '   -> ZoneType:',zone,zoneType!, Structured, Unstructured
  !      end if

  !      if(zoneType .ne. Unstructured)then
  !         print *,"Unstructured zone expected...stopping"
  !         stop
  !      end if

  !      ! Determine the number of vertices and volume elements in this 
  !      ! zone (and thus in the grid, because one zone is assumed). 
  !      call cg_zone_read_f(fileIndex,base,zone,zoneName,sizes,ierr)
  !      nVertices    = sizes(1)
  !      nVolElements = sizes(2)
  !      zones(zone)%nVertices = nVertices
  !      zones(zone)%nVolElements = nVolElements
  !      zones(zone)%zoneName = zoneName
  !      ! Allocate the memory for the grid points
  !      allocate(zones(zone)%points(nVertices,physDim),STAT=ierr)
  !      allocate(zones(zone)%points0(nVertices,physDim),STAT=ierr)
  !      allocate(zones(zone)%dx(nVertices,physDim),STAT=ierr)
  !      allocate(zones(zone)%isSurfaceNode(nVertices),STAT=ierr)
  !      zones(zone)%isSurfaceNode=.False.
  !      ! Read the x,y,z-coordinates.
  !      ! This assumes the coords are in double precision.
  !      ! Note that CGNS starts the numbering at 
  !      ! 1 even if C is used. 

  !      call cg_coord_read_f(fileIndex,base,zone,"CoordinateX",&
  !           & realDouble, 1,nVertices, zones(zone)%points(:,1)&
  !           &,ierr)

  !      call cg_coord_read_f(fileIndex,base,zone,"CoordinateY",&
  !           & realDouble, 1,nVertices, zones(zone)%points(:,2),&
  !           &ierr)

  !      call cg_coord_read_f(fileIndex,base,zone,"CoordinateZ",&
  !           & realDouble, 1,nVertices, zones(zone)%points(:,3),&
  !           &ierr)

  !      ! copy points to points0 to keep original points as a reference
  !      zones(zone)%points0 = zones(zone)%points
  !      ! Determine the number of sections for this zone. Note that 
  !      ! surface elements can be stored in a volume zone, but they 
  !      ! are NOT taken into account in the number obtained from 
  !      ! cg_zone_read. 

  !      call cg_nsections_f(fileIndex,base,zone,nSections,ierr)
  !      if (myID == 0) then
  !         print *, '   -> nSections',nSections
  !      end if
  !      ! zones(zone)%nSections = nSections

  !      ! ! allocate the logical to track the section type
  !      ! allocate(zones(zone)%isVolumeSection(nSections),STAT=ierr)

  !      call createSections(zone,nSections)

  !      ! loop over the number of sections and determine which
  !      ! sections are volume elements and which are surface elements

  !      zones(zone)%nVolSec = 0
  !      zones(zone)%nSurfSec = 0

  !      do sec=1,nSections
  !         ! Determine the element type and set the pointer for the 
  !         ! connectivity accordingly. 
  !         call cg_section_read_f(fileIndex,base,zone,sec,secName,&
  !              type,eBeg,eEnd,nBdry,parentFlag,ierr)

  !         select case (type)

  !         case (TETRA_4)
  !            zones(zone)%nVolSec = zones(zone)%nVolSec+1
  !            zones(zone)%isVolumeSection(sec)=.true.
  !         case (PYRA_5)
  !            zones(zone)%nVolSec = zones(zone)%nVolSec+1
  !            zones(zone)%isVolumeSection(sec)=.true.
  !         case (PENTA_6)
  !            zones(zone)%nVolSec = zones(zone)%nVolSec+1
  !            zones(zone)%isVolumeSection(sec)=.true.
  !         case (HEXA_8)
  !            zones(zone)%nVolSec = zones(zone)%nVolSec+1
  !            zones(zone)%isVolumeSection(sec)=.true.
  !         case (TRI_3)
  !            zones(zone)%nSurfSec = zones(zone)%nSurfSec+1
  !            zones(zone)%isVolumeSection(sec)=.false.
  !         case (QUAD_4)
  !            zones(zone)%nSurfSec =  zones(zone)%nSurfSec+1
  !            zones(zone)%isVolumeSection(sec)=.false.
  !         case  default
  !            if (myID == 0) then
  !               print *, "Unsupported element encountered....exiting", type
  !               stop
  !            end if
  !         end select

  !      end do

  !      allocate(zones(zone)%volumeSections(zones(zone)%nVolSec),STAT=ierr)
  !      allocate(zones(zone)%surfaceSections(zones(zone)%nSurfSec),STAT=ierr)

  !      ! Loop over the number of sections and read the element 
  !      ! connectivities. As CGNS starts the numbering at 1 the 
  !      ! for-loop starts at 1 as well. 

  !      volSecCounter = 1
  !      surfSecCounter = 1

  !      do sec=1,nSections
  !         ! Determine the element type and set the pointer for the 
  !         ! connectivity accordingly. 
  !         call cg_section_read_f(fileIndex,base,zone,sec,secName,&
  !              type,eBeg,eEnd,nBdry,parentFlag,ierr)
  !         nElem = eEnd-eBeg+1
  !         !print *,'type',sec,secName,type,eBeg,eEnd,nBdry,parentFlag
  !         !print *,'Types',HEXA_8,QUAD_4

  !         call cg_npe_f(type, nconn, ierr)

  !         if (zones(zone)%isVolumeSection(sec)) then
  !            ! this is a group of volume elements
  !            allocate(zones(zone)%volumeSections(volSecCounter)%elemPtr(nElem+1),STAT=ierr)
  !            allocate(zones(zone)%volumeSections(volSecCounter)%elemConn(nConn*nElem),STAT=ierr)
  !            zones(zone)%volumeSections(volSecCounter)%nElem = nElem
  !            zones(zone)%volumeSections(volSecCounter)%secName = secName
  !            zones(zone)%volumeSections(volSecCounter)%elemType = type
  !            !zones(zone)%volumeSections(volSecCounter)%nConn = nConn 
  !            zones(zone)%volumeSections(volSecCounter)%elemStart = eBeg
  !            zones(zone)%volumeSections(volSecCounter)%elemEnd = eEnd
  !         else
  !            ! this is a group of surface elements
  !            allocate(zones(zone)%surfaceSections(surfSecCounter)%elemPtr(nElem+1),STAT=ierr)
  !            allocate(zones(zone)%surfaceSections(surfSecCounter)%elemConn(nConn*nElem),STAT=ierr)

  !            ! allocate(zones(zone)%surfaceSections(surfSecCounter)%elements(nConn,nElem),STAT=ierr)
  !            allocate(zones(zone)%surfaceSections(surfSecCounter)%elemCenter(nElem,physDim),STAT=ierr)
  !            allocate(zones(zone)%surfaceSections(surfSecCounter)%elemArea(nElem,physDim),STAT=ierr)
  !            allocate(zones(zone)%surfaceSections(surfSecCounter)%elemNormal(nElem,physDim),STAT=ierr)
  !            allocate(zones(zone)%surfaceSections(surfSecCounter)%elemAreaMag(nElem),STAT=ierr)
  !            zones(zone)%surfaceSections(surfSecCounter)%nElem = nElem
  !            zones(zone)%surfaceSections(surfSecCounter)%secName = secName
  !            zones(zone)%surfaceSections(surfSecCounter)%elemType = type
  !            !zones(zone)%surfaceSections(surfSecCounter)%nConn = nConn 
  !            zones(zone)%surfaceSections(surfSecCounter)%elemStart = eBeg
  !            zones(zone)%surfaceSections(surfSecCounter)%elemEnd = eEnd        
  !         end if

  !         ! call cg_npe_f(type, nconn, ierr)
  !         allocate(conn(nConn,nElem),STAT=ierr)
  !         conn = 0
  !         ! allocate(parentData(nElem,nConn),STAT=ierr)
  !         ! Read the connectivity. Again, the node numbering of the 
  !         ! connectivities start at 1. If internally a starting
  !         ! index of 0 is used (typical for C-codes) 1 must be
  !         ! substracted from the connectivities read. 

  !         call cg_elements_read_f(fileIndex,base,zone,sec,conn,parentFlag,ierr)
  !         ! if (myID == 0) then
  !         !    print *, 'conn',shape(conn),conn(:,1)
  !         ! end if

  !         if (zones(zone)%isVolumeSection(sec)) then
  !            counter = 1
  !            do i = 1,nElem
  !               zones(zone)%volumeSections(volSecCounter)%elemPtr(i)=counter
  !               do j = 1,nConn
  !                  zones(zone)%volumeSections(volSecCounter)%elemConn(counter)=conn(j,i)
  !                  counter = counter +1
  !                  !zones(zone)%volumeSections(volSecCounter)%elements = conn
  !               end do
  !            end do
  !            zones(zone)%volumeSections(volSecCounter)%elemPtr(i)=counter
  !            volSecCounter = volSecCounter+1
  !         else
  !            !zones(zone)%surfaceSections(surfSecCounter)%elements = conn
  !            counter = 1
  !            do i = 1,nElem
  !               zones(zone)%surfaceSections(surfSecCounter)%elemPtr(i)=counter
  !               do j = 1,nConn
  !                  zones(zone)%surfaceSections(surfSecCounter)%elemConn(counter)=conn(j,i)
  !                  counter = counter +1
  !               end do
  !               !print *,'counter',i,counter,zones(zone)%surfaceSections(surfSecCounter)%elemPtr(i)
  !            end do
  !            !print *,'end',i,counter
  !            zones(zone)%surfaceSections(surfSecCounter)%elemPtr(i)=counter
  !            surfSecCounter = surfSecCounter + 1
  !         end if
  !         deallocate(conn,STAT=ierr)
  !      end do

  !      !Determine the number of boundary conditions for this zone. 
  !      call cg_nbocos_f(fileIndex,base,zone,nBocos,ierr)

  !      ! Loop over the number of boundary conditions.
  !      zones(zone)%nBocos = nBocos
  !      allocate(zones(zone)%BCInfo(nBocos),STAT=ierr)

  !      do boco=1,nBocos
  !         ! Read the info for this boundary condition. 

  !         call cg_boco_info_f(fileIndex, base, zone, boco, bocoName,&
  !              bocoType, ptsettype, nbcelem, normalIndex,&
  !              normListFlag,normDataType,nDataSet,ierr)
  !         ! print *,'boco info', boco, bocoName,&
  !         !      bocoType, ptsettype, nbcelem, normalIndex,&
  !         !      normListFlag,normDataType,nDataSet

  !         zones(zone)%BCInfo(boco)%nBCElem = nBCElem
  !         zones(zone)%BCInfo(boco)%BCName =  bocoName
  !         zones(zone)%BCInfo(boco)%BCType = bocoType
  !         allocate(zones(zone)%BCInfo(boco)%BCElements(nBCElem),STAT=ierr)

  !         ! Read the element ID's.
  !         call cg_boco_read_f(fileIndex,base,zone,boco,&
  !              zones(zone)%BCInfo(boco)%BCElements,NULL,ierr)

  !         ! Now that we have the boundary condition info read the family name
  !         ! to associate it with a specific section

  !         call cg_goto_f(fileIndex, base, ierr, "Zone_t", zone, "ZoneBC_t", 1, "BC_t", boco, "end")

  !         if (ierr == 0) then ! Node exits
  !            call cg_famname_read_f(famName, ierr)
  !            if (ierr == 0) then
  !               !search for matching section
  !               ! Loop over the number of sections and compare the family
  !               ! name to the section name. If they match, we have found
  !               ! the section that this BC belongs to.

  !               volSecCounter = 1
  !               surfSecCounter = 1
  !               do sec=1,nSections
  !                  if (.not. zones(zone)%isVolumeSection(sec)) then
  !                     ! we have a surface section
  !                     zones(zone)%BCInfo(boco)%famName = famName
  !                     !print *,'names',famName,zones(zone)%surfaceSections(surfSecCounter)%secName
  !                     if(zones(zone)%surfaceSections(surfSecCounter)%secName==&
  !                          famName)then
  !                        !print *,'names Match',famName,zones(zone)%surfaceSections(surfSecCounter)%secName
  !                        ! this is the correct section
  !                        ! Set BCType and Family Name
  !                        zones(zone)%surfaceSections(surfSecCounter)%BCFamily = famName
  !                        zones(zone)%surfaceSections(surfSecCounter)%BCType = bocoType
  !                     ! else
  !                     !    ! this is not the correct section
  !                     end if
  !                     surfSecCounter = surfSecCounter + 1
  !                  end if
  !               end do

  !               !print *,'BCs',zone,index,zones(zone)%surfaceSections(index)%BCFamily,zones(zone)%surfaceSections(index)%BCType
  !            ! else
  !            !    !family not present set logical
  !            end if
  !         end if

  !      end do

  !   end do


  !  ! We can also generate the wall family list
  !   nwallFamilies = 0
  !   familyList(:) = ''

  !   do zone=1, nzones
  !      surfSecCounter = 1
  !      do sec=1,nSections
  !         if (.not. zones(zone)%isVolumeSection(sec)) then
  !            bocoType = zones(zone)%surfaceSections(surfSecCounter)%BCType
  !            FamName = zones(zone)%surfaceSections(surfSecCounter)%BCFamily
  !            zones(zone)%surfaceSections(surfSecCounter)%isWallBC=.False.
  !            zones(zone)%surfaceSections(surfSecCounter)%isBoundaryBC=.False.
  !            zones(zone)%surfaceSections(surfSecCounter)%isSymmBC=.False.
  ! !           print *,'bocotype',bocoType,BCTypeName(bocoType)
  !            if (BCTypeName(bocoType) == 'BCWallViscous' .or. &
  !                 BCTypeName(bocoType) == 'BCWallInviscid' .or. &
  !                 BCTypeName(bocoType) == 'BCWall' .or. &
  !                 BCTypeName(bocoType) == 'BCWallViscousHeatFlux' .or. &
  !                 BCTypeName(bocoType) == 'BCWallViscousIsothermal') then
  !               !print *,'is Wall',zone,surfSecCounter
  !               zones(zone)%surfaceSections(surfSecCounter)%isWallBC=.True.
  !               !this is a wall family, add to the family list
  !               call checkInFamilyList(familyList, famName, nwallFamilies, famID)
  !               if (famID == 0) then
  !                  nwallFamilies = nwallFamilies + 1
  !                  familyList(nwallFamilies) = famName
  !               end if
  !            elseif(BCTypeName(bocotype) == 'BCFarfield')  then
  !               !print *,'is farfield',zone,surfSecCounter
  !               zones(zone)%surfaceSections(surfSecCounter)%isBoundaryBC=.True.
  !            elseif(BCTypeName(bocotype) == 'BCSymmetryPlane') then
  !               !print *,' is symm',zone,surfSecCounter
  !               zones(zone)%surfaceSections(surfSecCounter)%isSymmBC = .True.
  !               hasSymmetry = .True.
  !            else
  !               print *,'Unrecongnized boundary Type... exiting'
  !               stop
  !            end if
  !            surfSecCounter = surfSecCounter + 1
  !         end if
  !      end do
  !   end do

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


