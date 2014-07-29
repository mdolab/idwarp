! ====================================================================
! File: writeUnstructuredCGNS.f90
! Author: C.A.(Sandy) Mader
! Date Started: July 10, 2014
! Date Modified:

subroutine writeUnstructuredCGNSFile(cgns_file, comm)

  use precision
  use communication
  use gridData
  use cgnsData
  implicit none
  !include 'cgnslib_f.h'

  ! Input Arguments
  character*(*) :: cgns_file
  integer(kind=intType) :: comm

  integer(kind=intType):: ierr

  character*32 :: zoneName,secName,bocoName,famName

  ! ! CGNS Variabls
  integer(kind=intType) :: fileIndex, base
  integer(kind=intType) :: nbocos, boco,bcIndex
  integer(kind=intType) :: nVertices,nVolElements
  integer(kind=intType) :: zone,sec,type
  integer(kind=intType) :: coordIndex
  integer(kind=intType) :: nSections,nConn,nElem
  integer(kind=intType) :: eBeg,eEnd
  integer(kind=intType) :: bocoType,nbcelem, ipnts
  integer(kind=intType), dimension(:,:), allocatable::conn

  integer(kind=intType):: surfSecCounter,volSecCounter

  ! ---------------------------------------
  !           Open CGNS File
  ! ---------------------------------------
  if (myID == 0) then
     print *, ' -> Writing CGNS File: ', cgns_file
  end if

  ! Open the CGNS File
  call cg_open_f(cgns_file, CG_MODE_WRITE, fileIndex, ierr)
  if (ierr .eq. CG_ERROR) call cg_error_exit_f

  ! Now write the base entry in the file
  call cg_base_write_f(fileIndex,baseName,celldim,physdim,&
       base,ierr)

  ! loop over the zones, create and write
  do zone = 1,nZones
     ! Get the name and number of points and elements that define the zone
     zoneName = gridDoms(zone)%zoneName
     sizes(1) = gridDoms(zone)%nVertices
     sizes(2) = gridDoms(zone)%nVolElements

     ! Create the zone
     call cg_zone_write_f(fileIndex,base,zoneName,sizes,&
          Unstructured,zone,ierr)

     ! Write the x,y,z-coordinates.
     ! This assumes the coords are in double precision.
     ! Note that CGNS starts the numbering at 
     ! 1 even if C is used. 

     call cg_coord_write_f(fileIndex,base,zone,realDouble,&
          "CoordinateX",gridDoms(zone)%points(:,1),coordIndex,&
          ierr)

     call cg_coord_write_f(fileIndex,base,zone,realDouble,&
          "CoordinateY",gridDoms(zone)%points(:,2),coordIndex,&
          ierr)

     call cg_coord_write_f(fileIndex,base,zone, realDouble,&
          "CoordinateZ",gridDoms(zone)%points(:,3),coordIndex,&
          ierr)
     nSections = gridDoms(zone)%nSections
     volSecCounter = 1
     surfSecCounter = 1
     do sec=1,nSections
        ! Get the section information for the current section
        if (gridDoms(zone)%isVolumeSection(sec)) then
           nElem = gridDoms(zone)%volumeSections(volSecCounter)%nElem
           secName = gridDoms(zone)%volumeSections(volSecCounter)%secName
           type = gridDoms(zone)%volumeSections(volSecCounter)%elemType
           nConn = gridDoms(zone)%volumeSections(volSecCounter)%nConn
           eBeg = gridDoms(zone)%volumeSections(volSecCounter)%elemStart
           eEnd = gridDoms(zone)%volumeSections(volSecCounter)%elemEnd
           allocate(conn(nConn,nElem),STAT=ierr)
           conn = gridDoms(zone)%volumeSections(volSecCounter)%elements
           volSecCounter = volSecCounter+1
        else
           nElem = gridDoms(zone)%surfaceSections(surfSecCounter)%nElem
           secName =  gridDoms(zone)%surfaceSections(surfSecCounter)%secName
           type = gridDoms(zone)%surfaceSections(surfSecCounter)%elemType
           nConn = gridDoms(zone)%surfaceSections(surfSecCounter)%nConn
           eBeg = gridDoms(zone)%surfaceSections(surfSecCounter)%elemStart
           eEnd = gridDoms(zone)%surfaceSections(surfSecCounter)%elemEnd
           allocate(conn(nConn,nElem),STAT=ierr)
           conn = gridDoms(zone)%surfaceSections(surfSecCounter)%elements
           surfSecCounter = surfSecCounter + 1
        end if

        !write element connectivity (user can give any name)
        call cg_section_write_f(fileIndex,base,zone,secName,type,&
             eBeg,eEnd,0,conn,sec,ierr)
        deallocate(conn,STAT=ierr)

     end do

     ! Write the BCs for the mesh
     nBocos = gridDoms(zone)%nBocos

     do boco=1,nBocos

        ! retrieve the boundary data
        nBCElem = gridDoms(zone)%BCInfo(boco)%nBCElem
        bocoName = gridDoms(zone)%BCInfo(boco)%BCName
        bocoType = gridDoms(zone)%BCInfo(boco)%BCType

        ! Write each of the boundary conditions. 
        call cg_boco_write_f(fileIndex,base,zone,bocoName,bocoType,&
             ElementList,nBCElem,gridDoms(zone)%BCInfo(boco)%BCElements,&
             bcIndex,ierr)

        ! ! now write the family Name
        call cg_goto_f(fileIndex, base, ierr, "Zone_t", zone, "ZoneBC_t", 1, "BC_t", boco, "end")
        
        if (ierr == 0) then ! Node exits
           call cg_famname_write_f(gridDoms(zone)%BCInfo(boco)%famName, ierr)
        end if

     end do

  end do
!  close CGNS file
  call cg_close_f(fileIndex,ierr)

end subroutine writeUnstructuredCGNSFile

