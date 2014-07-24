! ====================================================================
! File: readUnstructuredCGNS.f90
! Author: C.A.(Sandy) Mader
! Date Started: July 10, 2014
! Date Modified:

subroutine readUnstructuredCGNSFile(cgns_file, comm)

  use precision
  use communication
  use gridData
  implicit none
  include 'cgnslib_f.h'

  ! Input Arguments
  character*(*) :: cgns_file
  integer(kind=intType) :: comm

  integer(kind=intType):: ierr
  ! integer(kind=intType):: index_base, index_zone
  ! integer(kind=intType), dimension(3,1):: isize
  ! integer(kind=intType), dimension(3):: irmin, irmax

  character*32 :: zoneName
  character*32 :: baseName,name,secName,bocoName,famName
  ! CGNS Variabls
  integer(kind=intType) :: fileIndex, base, nbases
  integer(kind=intType) :: nbocos, boco, index
  integer(kind=intType) :: nVertices,nVolElements
  integer(kind=intType) :: zone,zoneType,dataType,sec,type
  integer(kind=intType) :: nCoords,nSections,nConn,nElem
  integer(kind=intType) :: eBeg,eEnd,nBdry, parentFlag
  integer(kind=intType) :: bocoType,ptsettype,nbcelem
  integer(kind=intType) :: normalIndex,normListFlag
  integer(kind=intType) :: normDataType,nDataSet
  !integer(kind=8), dimension(3)::sizes
  integer(kind=intType), dimension(3)::sizes
  integer(kind=intType), dimension(:,:), allocatable::conn

  integer(kind=intType):: surfSecCounter,volSecCounter, famID
  !integer(kind=intType), dimension(:), allocatable:: BCElements
  !integer, dimension(:), allocatable::conn,parentData

  ! ---------------------------------------
  !            Setup MPI Data first
  ! ---------------------------------------

  warp_comm_world = comm
  warp_comm_self  = mpi_comm_self
  call MPI_Comm_size(warp_comm_world, nProc, ierr)
  call MPI_Comm_rank(warp_comm_world, myid , ierr)

  ! ---------------------------------------
  !           Open CGNS File
  ! ---------------------------------------
  if (myID == 0) then
     print *, ' -> Reading CGNS File: ', cgns_file
  end if

  ! if (CG_BUILD_64BIT .ne. 0) then
  !    print *,'will not work in 64-bit mode',CG_BUILD_64BIT 
  !    stop
  ! endif

  call cg_open_f(cgns_file, CG_MODE_READ, fileIndex, ierr)
  if (ierr .eq. CG_ERROR) call cg_error_exit_f

  call cg_nbases_f(fileIndex, nbases, ierr)
  if (ierr .eq. CG_ERROR) call cg_error_exit_f

  if (nbases .gt. 1) then
     print *, ' ** Warning: pyWarp only reads the first base in a cgns file'
  end if

  base = 1_intType

  !       *** base attribute:  GOTO base node
  call cg_goto_f(fileIndex, base, ierr, 'end')
  if (ierr .eq. CG_ERROR) call cg_error_exit_f

  ! Check the cell and physical dimensions of the bases. 
  ! Both should be 3. 
  call cg_base_read_f(fileIndex,base,baseName,cellDim,physDim,ierr)

  if (ierr .eq. ERROR) call cg_error_exit_f
  if (myID == 0) then
     print *, '   -> BaseName and Dimensions:', baseName,cellDim, physDim
  end if

  if( cellDim .ne.3 .or. physDim .ne. 3) then
     print *,' ** Warning:  pyWarp only accepts 3-d data'
  end if

  ! read the number of zones in the file
  call cg_nzones_f(fileIndex, base, nzones, ierr)
  if (ierr .eq. ERROR) call cg_error_exit_f
  if (myID == 0) then
     print *, '   -> Number of Zones:', nzones
  end if

  ! Allocate the GridDoms Data structure
  allocate(gridDoms(nZones),STAT=ierr)

  ! loop over the zones and read
  do zone = 1,nZones

     ! Check the zone type. This should be Unstructured. 
     call cg_zone_type_f(fileIndex,base, zone, zoneType,ierr)
     if (ierr .eq. ERROR) call cg_error_exit_f
     if (myID == 0) then
        print *, '   -> ZoneType:',zone,zoneType!, Structured, Unstructured
     end if

     if(zoneType .ne. Unstructured)then
        print *,"Unstructured zone expected...stopping"
        stop
     end if

     ! Determine the number of vertices and volume elements in this 
     ! zone (and thus in the grid, because one zone is assumed). 
     call cg_zone_read_f(fileIndex,base,zone,zoneName,sizes,ierr)
     nVertices    = sizes(1)
     nVolElements = sizes(2)
     gridDoms(zone)%nVertices = nVertices

     ! Allocate the memory for the grid points
     allocate(gridDoms(zone)%points(nVertices,physDim),STAT=ierr)
     allocate(gridDoms(zone)%points0(nVertices,physDim),STAT=ierr)
     allocate(gridDoms(zone)%dx(nVertices,physDim),STAT=ierr)
     allocate(gridDoms(zone)%isSurfaceNode(nVertices),STAT=ierr)
     gridDoms(zone)%isSurfaceNode=.False.
     ! Read the x,y,z-coordinates.
     ! This assumes the coords are in double precision.
     ! Note that CGNS starts the numbering at 
     ! 1 even if C is used. 

     call cg_coord_read_f(fileIndex,base,zone,"CoordinateX",&
          & realDouble, 1,nVertices, gridDoms(zone)%points(:,1)&
          &,ierr)

     call cg_coord_read_f(fileIndex,base,zone,"CoordinateY",&
          & realDouble, 1,nVertices, gridDoms(zone)%points(:,2),&
          &ierr)

     call cg_coord_read_f(fileIndex,base,zone,"CoordinateZ",&
          & realDouble, 1,nVertices, gridDoms(zone)%points(:,3),&
          &ierr)

     ! copy points to points0 to keep original points as a reference
     gridDoms(zone)%points0 = gridDoms(zone)%points
     ! Determine the number of sections for this zone. Note that 
     ! surface elements can be stored in a volume zone, but they 
     ! are NOT taken into account in the number obtained from 
     ! cg_zone_read. 

     call cg_nsections_f(fileIndex,base,zone,nSections,ierr)
     if (myID == 0) then
        print *, '   -> nSections',nSections
     end if
     gridDoms(zone)%nSections = nSections

     ! allocate the logical to track the section type
     allocate(gridDoms(zone)%isVolumeSection(nSections),STAT=ierr)

     ! loop over the number of sections and determine which
     ! sections are volume elements and which are surface elements

     gridDoms(zone)%nVolSec = 0
     gridDoms(zone)%nSurfSec = 0

     do sec=1,nSections
        ! Determine the element type and set the pointer for the 
        ! connectivity accordingly. 
        call cg_section_read_f(fileIndex,base,zone,sec,secName,&
             type,eBeg,eEnd,nBdry,parentFlag,ierr)

        select case (type)

        case (TETRA_4)
           gridDoms(zone)%nVolSec = gridDoms(zone)%nVolSec+1
           gridDoms(zone)%isVolumeSection(sec)=.true.
        case (PYRA_5)
           gridDoms(zone)%nVolSec = gridDoms(zone)%nVolSec+1
           gridDoms(zone)%isVolumeSection(sec)=.true.
        case (PENTA_6)
           gridDoms(zone)%nVolSec = gridDoms(zone)%nVolSec+1
           gridDoms(zone)%isVolumeSection(sec)=.true.
        case (HEXA_8)
           gridDoms(zone)%nVolSec = gridDoms(zone)%nVolSec+1
           gridDoms(zone)%isVolumeSection(sec)=.true.
        case (TRI_3)
           gridDoms(zone)%nSurfSec = gridDoms(zone)%nSurfSec+1
           gridDoms(zone)%isVolumeSection(sec)=.false.
        case (QUAD_4)
           gridDoms(zone)%nSurfSec =  gridDoms(zone)%nSurfSec+1
           gridDoms(zone)%isVolumeSection(sec)=.false.
        case  default
           if (myID == 0) then
              print *, "Unsupported element encountered....exiting", type
              stop
           end if
        end select

     end do

     allocate(gridDoms(zone)%volumeSections(gridDoms(zone)%nVolSec),STAT=ierr)
     allocate(gridDoms(zone)%surfaceSections(gridDoms(zone)%nSurfSec),STAT=ierr)

     ! Loop over the number of sections and read the element 
     ! connectivities. As CGNS starts the numbering at 1 the 
     ! for-loop starts at 1 as well. 

     volSecCounter = 1
     surfSecCounter = 1

     do sec=1,nSections
        ! Determine the element type and set the pointer for the 
        ! connectivity accordingly. 
        call cg_section_read_f(fileIndex,base,zone,sec,secName,&
             type,eBeg,eEnd,nBdry,parentFlag,ierr)
        nElem = eEnd-eBeg+1
        !print *,'type',sec,secName,type,eBeg,eEnd,nBdry,parentFlag
        !print *,'Types',HEXA_8,QUAD_4

        call cg_npe_f(type, nconn, ierr)

        if (gridDoms(zone)%isVolumeSection(sec)) then
           ! this is a group of volume elements
           allocate(gridDoms(zone)%volumeSections(volSecCounter)%elements(nConn,nElem),STAT=ierr)
           gridDoms(zone)%volumeSections(volSecCounter)%nElem = nElem
           gridDoms(zone)%volumeSections(volSecCounter)%secName = secName
           gridDoms(zone)%volumeSections(volSecCounter)%elemType = type
           gridDoms(zone)%volumeSections(volSecCounter)%nConn = nConn 
           gridDoms(zone)%volumeSections(volSecCounter)%elemStart = eBeg
           gridDoms(zone)%volumeSections(volSecCounter)%elemEnd = eEnd
        else
           ! this is a group of surface elements
           allocate(gridDoms(zone)%surfaceSections(surfSecCounter)%elements(nConn,nElem),STAT=ierr)
           allocate(gridDoms(zone)%surfaceSections(surfSecCounter)%elemCenter(nElem,physDim),STAT=ierr)
           allocate(gridDoms(zone)%surfaceSections(surfSecCounter)%elemArea(nElem,physDim),STAT=ierr)
           allocate(gridDoms(zone)%surfaceSections(surfSecCounter)%elemNormal(nElem,physDim),STAT=ierr)
           allocate(gridDoms(zone)%surfaceSections(surfSecCounter)%elemAreaMag(nElem),STAT=ierr)
           gridDoms(zone)%surfaceSections(surfSecCounter)%nElem = nElem
           gridDoms(zone)%surfaceSections(surfSecCounter)%secName = secName
           gridDoms(zone)%surfaceSections(surfSecCounter)%elemType = type
           gridDoms(zone)%surfaceSections(surfSecCounter)%nConn = nConn 
           gridDoms(zone)%surfaceSections(surfSecCounter)%elemStart = eBeg
           gridDoms(zone)%surfaceSections(surfSecCounter)%elemEnd = eEnd        
        end if

        ! call cg_npe_f(type, nconn, ierr)
        allocate(conn(nConn,nElem),STAT=ierr)
        conn = 0
        ! allocate(parentData(nElem,nConn),STAT=ierr)
        ! Read the connectivity. Again, the node numbering of the 
        ! connectivities start at 1. If internally a starting
        ! index of 0 is used (typical for C-codes) 1 must be
        ! substracted from the connectivities read. 

        call cg_elements_read_f(fileIndex,base,zone,sec,conn,parentFlag,ierr)
        ! if (myID == 0) then
        !    print *, 'conn',shape(conn),conn(:,1)
        ! end if

        if (gridDoms(zone)%isVolumeSection(sec)) then
           gridDoms(zone)%volumeSections(volSecCounter)%elements = conn
           volSecCounter = volSecCounter+1
        else
           gridDoms(zone)%surfaceSections(surfSecCounter)%elements = conn
           surfSecCounter = surfSecCounter + 1
        end if
        deallocate(conn,STAT=ierr)
     end do

     ! ignore BCs for now. All the info we need is in the sections

     !Determine the number of boundary conditions for this zone. 
     call cg_nbocos_f(fileIndex,base,zone,nBocos,ierr)

     ! Loop over the number of boundary conditions.

     do boco=1,nBocos
        ! Read the info for this boundary condition. 

        call cg_boco_info_f(fileIndex, base, zone, boco, bocoName,&
             bocoType, ptsettype, nbcelem, normalIndex,&
             normListFlag,normDataType,nDataSet,ierr)
        ! print *,'boco info', boco, bocoName,&
        !      bocoType, ptsettype, nbcelem, normalIndex,&
        !      normListFlag,normDataType,nDataSet

        ! Now that we have the boundary condition info read the family name
        ! to associate it with a specific section

        call cg_goto_f(fileIndex, base, ierr, "Zone_t", zone, "ZoneBC_t", 1, "BC_t", boco, "end")
        
        if (ierr == 0) then ! Node exits
           call cg_famname_read_f(famName, ierr)
           if (ierr == 0) then
              !search for matching section
              ! Loop over the number of sections and compare the family
              ! name to the section name. If they match, we have found
              ! the section that this BC belongs to.
              
              volSecCounter = 1
              surfSecCounter = 1
              do sec=1,nSections
                 if (.not. gridDoms(zone)%isVolumeSection(sec)) then
                    ! we have a surface section
                    !print *,'names',famName,gridDoms(zone)%surfaceSections(surfSecCounter)%secName
                    if(gridDoms(zone)%surfaceSections(surfSecCounter)%secName==&
                         famName)then
                       !print *,'names Match',famName,gridDoms(zone)%surfaceSections(surfSecCounter)%secName
                       ! this is the correct section
                       ! Set BCType and Family Name
                       gridDoms(zone)%surfaceSections(surfSecCounter)%BCFamily = famName
                       gridDoms(zone)%surfaceSections(surfSecCounter)%BCType = bocoType
                    ! else
                    !    ! this is not the correct section
                    end if
                    surfSecCounter = surfSecCounter + 1
                 end if
              end do
            
              !print *,'BCs',zone,index,gridDoms(zone)%surfaceSections(index)%BCFamily,gridDoms(zone)%surfaceSections(index)%BCType
           ! else
           !    !family not present set logical
           end if
        end if
        ! !         allocate(BCElements(nbcelem),STAT=ierr)
        ! !         ! Read the element ID's.
        ! !         call cg_boco_read_f(fileIndex,base,zone,boco,&
        ! !              BCElements,NULL,ierr)
        ! ! !        print *,'bcelem',shape(BCElements)
        ! !         deallocate(BCElements,STAT=ierr)
     end do

  end do

 ! We can also generate the wall family list
  nwallFamilies = 0
  familyList(:) = ''
  
  do zone=1, nzones
     surfSecCounter = 1
     do sec=1,nSections
        if (.not. gridDoms(zone)%isVolumeSection(sec)) then
           bocoType = gridDoms(zone)%surfaceSections(surfSecCounter)%BCType
           FamName = gridDoms(zone)%surfaceSections(surfSecCounter)%BCFamily
           gridDoms(zone)%surfaceSections(surfSecCounter)%isWallBC=.False.
           gridDoms(zone)%surfaceSections(surfSecCounter)%isBoundaryBC=.False.
           gridDoms(zone)%surfaceSections(surfSecCounter)%isSymmBC=.False.
           print *,'bocotype',bocoType,BCTypeName(bocoType)
           if (BCTypeName(bocoType) == 'BCWallViscous' .or. &
                BCTypeName(bocoType) == 'BCWallInviscid' .or. &
                BCTypeName(bocoType) == 'BCWall' .or. &
                BCTypeName(bocoType) == 'BCWallViscousHeatFlux' .or. &
                BCTypeName(bocoType) == 'BCWallViscousIsothermal') then
              print *,'is Wall',zone,surfSecCounter
              gridDoms(zone)%surfaceSections(surfSecCounter)%isWallBC=.True.
              !this is a wall family, add to the family list
              call checkInFamilyList(familyList, famName, nwallFamilies, famID)
              if (famID == 0) then
                 nwallFamilies = nwallFamilies + 1
                 familyList(nwallFamilies) = famName
              end if
           elseif(BCTypeName(bocotype) == 'BCFarfield')  then
              print *,'is farfield',zone,surfSecCounter
              gridDoms(zone)%surfaceSections(surfSecCounter)%isBoundaryBC=.True.
           elseif(BCTypeName(bocotype) == 'BCSymmetryPlane') then
              print *,' is symm',zone,surfSecCounter
              gridDoms(zone)%surfaceSections(surfSecCounter)%isSymmBC = .True.
           else
              print *,'Unrecongnized boundary Type... exiting'
              stop
           end if
           surfSecCounter = surfSecCounter + 1
        end if
     end do
  end do

end subroutine readUnstructuredCGNSFile


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

subroutine deallocateCGNSData()

  use precision
  use communication
  use gridData
  implicit none

  integer(kind=intType)::ierr
  integer(kind=intType)::zone,sec
  integer(kind=intType)::surfSecCounter,volSecCounter
  
  ! loop over the zones and deallocate Memory
  do zone = 1,nZones
     volSecCounter = 1
     surfSecCounter = 1
     do sec=1,gridDoms(zone)%nSections
        
        if (gridDoms(zone)%isVolumeSection(sec)) then
           deallocate(gridDoms(zone)%volumeSections(volSecCounter)%elements,STAT=ierr)
            volSecCounter = volSecCounter + 1
        else
           deallocate(gridDoms(zone)%surfaceSections(surfSecCounter)%elements,STAT=ierr)
           deallocate(gridDoms(zone)%surfaceSections(surfSecCounter)%elemCenter,STAT=ierr)
           deallocate(gridDoms(zone)%surfaceSections(surfSecCounter)%elemArea,STAT=ierr)
           deallocate(gridDoms(zone)%surfaceSections(surfSecCounter)%elemNormal,STAT=ierr)
           deallocate(gridDoms(zone)%surfaceSections(surfSecCounter)%elemAreaMag,STAT=ierr)
           surfSecCounter = surfSecCounter + 1
        end if
     end do
     deallocate(gridDoms(zone)%volumeSections,STAT=ierr)
     deallocate(gridDoms(zone)%surfaceSections,STAT=ierr)
     deallocate(gridDoms(zone)%isVolumeSection,STAT=ierr)
     deallocate(gridDoms(zone)%dx,STAT=ierr)
     deallocate(gridDoms(zone)%points0,STAT=ierr)
     deallocate(gridDoms(zone)%points,STAT=ierr)
  end do

  deallocate(gridDoms,STAT=ierr)

end subroutine deallocateCGNSData
