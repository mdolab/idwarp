! **************************************************
! Routines for dealing with the structured CGNS mesh
! **************************************************

subroutine getNPatchesStructured(nPatches)

  ! Return the number of wall patches we have on this proc
  use precision
  use constants
  use structuredCGNSGrid
 
  implicit none
  include 'cgnslib_f.h'

  ! Output Variables
  integer(kind=intType),intent(out) :: nPatches

  ! Working Variables
  integer(kind=intType) :: zone, boco
  integer(kind=intType) :: nZones, nBCs
  integer(kind=intType) :: bocotype

  nPatches = 0_intType
  nZones = size(blocks) 

  do zone = 1,nZones
     nBCs = size(blocks(zone)%bocos) 
     do boco=1,nBCs
        bocoType = blocks(zone)%bocos(boco)%bocoType
        if (BCTypeName(bocoType) == 'BCWallViscous' .or. &
             BCTypeName(bocoType) == 'BCWallInviscid' .or. &
             BCTypeName(bocoType) == 'BCWall' .or. &
             BCTypeName(bocoType) == 'BCWallViscousHeatFlux' .or. &
             BCTypeName(bocoType) == 'BCWallViscousIsothermal') then
           ! this is a wall surface, count the patch
           nPatches = nPatches + 1
        end if
     end do
  end do
end subroutine getNPatchesStructured

subroutine getPatchNameStructured(iPatch, patchName) 
  ! Get the name of patch 'iPatch' one at a time since f2py doesn't
  !  like array's of strings.
  use precision
  use constants
  use structuredCGNSGrid

  implicit none
  include 'cgnslib_f.h'

  ! Input Variables
  integer(kind=intType),intent(in) :: iPatch

  ! Output Variables 
  character*32, intent(out) :: patchName

  ! Working
  integer(kind=intType) :: patchCount, zone, boco, bocoType

  patchCount = 0
  do zone = 1, size(blocks) 
     do boco=1, size(blocks(zone)%bocos)
        bocoType = blocks(zone)%bocos(boco)%bocoType
        if (BCTypeName(bocoType) == 'BCWallViscous' .or. &
             BCTypeName(bocoType) == 'BCWallInviscid' .or. &
             BCTypeName(bocoType) == 'BCWall' .or. &
             BCTypeName(bocoType) == 'BCWallViscousHeatFlux' .or. &
             BCTypeName(bocoType) == 'BCWallViscousIsothermal') then
           ! this is a wall surface, count the patch
           if (patchCount == iPatch) then
              patchName =  blocks(zone)%bocos(boco)%family
           end if
           patchCount = patchCount + 1
        end if
     end do
  end do
end subroutine getPatchNameStructured

subroutine getPatchSizeStructured(iPatch, patchSize)
  ! Get the sizes of patch 'iPatch'

  use precision
  use constants 
  use structuredCGNSGrid

  implicit none
  include 'cgnslib_f.h'
  
  ! Input Variables
  integer(kind=intType),intent(in) :: iPatch

  ! Output Variables 
  integer(kind=intType), intent(out) :: patchSize(2)

  ! Working
  integer(kind=intType) :: patchCount, zone, boco, bocoType
  integer(kind=intType) :: ptRange(3,2)

  patchCount = 0

  ! loop over the zones and read
  do zone = 1, size(blocks) 
     do boco=1, size(blocks(zone)%bocos) 
        bocoType = blocks(zone)%bocos(boco)%bocoType
        if (BCTypeName(bocoType) == 'BCWallViscous' .or. &
             BCTypeName(bocoType) == 'BCWallInviscid' .or. &
             BCTypeName(bocoType) == 'BCWall' .or. &
             BCTypeName(bocoType) == 'BCWallViscousHeatFlux' .or. &
             BCTypeName(bocoType) == 'BCWallViscousIsothermal') then
           
           ! this is a wall surface, count the patch
           if (patchCount == iPatch) then
              ptRange = blocks(zone)%bocos(boco)%ptRange
           
              if (ptRange(3,1)==ptRange(3,2)) then !Zmin/Zmax
                 patchSize(1) = blocks(zone)%bocos(boco)%ptRange(1,2)
                 patchSize(2) = blocks(zone)%bocos(boco)%ptRange(2,2)
              else if (ptRange(1,1)==ptRange(1,2)) then ! iMin/iMax
                 patchSize(1) = blocks(zone)%bocos(boco)%ptRange(2,2)
                 patchSize(2) = blocks(zone)%bocos(boco)%ptRange(3,2)
              else ! jMin, jMax
                 patchSize(1) = blocks(zone)%bocos(boco)%ptRange(1,2)
                 patchSize(2) = blocks(zone)%bocos(boco)%ptRange(3,2)
              end if
           end if
           patchCount = patchCount + 1
        end if
     end do
  end do
end subroutine getPatchSizeStructured

subroutine getNPatchesUnstructured(nPatches)

  ! Return the number of wall patches we have on this proc
  use precision
  use constants
  use unStructuredCGNSGrid
 
  implicit none
  include 'cgnslib_f.h'

  ! Output Variables
  integer(kind=intType),intent(out) :: nPatches

  ! Working Variables
  integer(kind=intType) :: iZone, sec
  integer(kind=intType) :: nSections
  integer(kind=intType) :: bocotype
  integer(kind=intType) :: surfSecCounter

  nPatches = 0_intType

  !nZones = size(zones) 

  ! loop over the zones and read
  do iZone = 1,nZones
     nSections = zones(iZone)%nSections
     surfSecCounter = 1
     do sec=1,nSections
        if (.not. zones(iZone)%isVolumeSection(sec)) then
           bocoType = zones(iZone)%surfaceSections(surfSecCounter)%BCType
           if (BCTypeName(bocoType) == 'BCWallViscous' .or. &
                BCTypeName(bocoType) == 'BCWallInviscid' .or. &
                BCTypeName(bocoType) == 'BCWall' .or. &
                BCTypeName(bocoType) == 'BCWallViscousHeatFlux' .or. &
                BCTypeName(bocoType) == 'BCWallViscousIsothermal') then
              ! this is a wall surface, count the patch
              nPatches = nPatches + 1
           end if
           surfSecCounter = surfSecCounter + 1
        end if
     end do
  end do
end subroutine getNPatchesUnstructured

! **************************************************
! Routines for dealing with the UNstructured CGNS mesh
! **************************************************

subroutine getPatchNameUnstructured(iPatch, patchName) 
  
  ! Get names one at a time since f2py can't do arrays of strings
  use precision
  use constants
  use unStructuredCGNSGrid

  implicit none
  include 'cgnslib_f.h'

  ! Input Variables
  integer(kind=intType),intent(in) :: iPatch

  ! Output Variables 
  character*32, intent(out) :: patchName

  ! Working
  integer(kind=intType) :: patchCount
  integer(kind=intType) :: iZone, sec
  integer(kind=intType) :: nSections
  integer(kind=intType) :: bocotype
  integer(kind=intType) :: surfSecCounter

  patchCount = 0
  ! loop over the zones and read
  do iZone = 1,nZones
     nSections = zones(iZone)%nSections
     surfSecCounter = 1       
     do sec=1,nSections
        if (.not. zones(iZone)%isVolumeSection(sec)) then
           bocoType = zones(iZone)%surfaceSections(surfSecCounter)%BCType
           if (BCTypeName(bocoType) == 'BCWallViscous' .or. &
                BCTypeName(bocoType) == 'BCWallInviscid' .or. &
                BCTypeName(bocoType) == 'BCWall' .or. &
                BCTypeName(bocoType) == 'BCWallViscousHeatFlux' .or. &
                BCTypeName(bocoType) == 'BCWallViscousIsothermal') then
              ! this is a wall surface, count the patch
              if (patchCount == iPatch) then
                 patchName =  zones(iZone)%surfaceSections(surfSecCounter)%BCFamily
              end if
              patchCount = patchCount + 1
           end if
           surfSecCounter = surfSecCounter + 1
        end if
     end do
  end do

end subroutine getPatchNameUnstructured


subroutine getPatchSizesUnstructured(iPatch, nPts, nCells) 
  
  ! Get names one at a time since f2py can't do arrays of strings
  use precision
  use constants
  use unStructuredCGNSGrid

  implicit none
  include 'cgnslib_f.h'

  ! Input Variables
  integer(kind=intType),intent(in) :: iPatch

  ! Output Variables 
  integer(kind=intType),intent(out) :: nPts,nCells


  ! Working
  integer(kind=intType) :: patchCount
  integer(kind=intType) :: iZone, sec
  integer(kind=intType) :: nSections
  integer(kind=intType) :: bocotype
  integer(kind=intType) :: surfSecCounter
  integer(kind=intType) :: nConn


  patchCount = 0
  ! loop over the zones and read
  do iZone = 1,nZones
     nSections = zones(iZone)%nSections
     surfSecCounter = 1       
     do sec=1,nSections
        if (.not. zones(iZone)%isVolumeSection(sec)) then
           bocoType = zones(iZone)%surfaceSections(surfSecCounter)%BCType
           if (BCTypeName(bocoType) == 'BCWallViscous' .or. &
                BCTypeName(bocoType) == 'BCWallInviscid' .or. &
                BCTypeName(bocoType) == 'BCWall' .or. &
                BCTypeName(bocoType) == 'BCWallViscousHeatFlux' .or. &
                BCTypeName(bocoType) == 'BCWallViscousIsothermal') then
              ! this is a wall surface, count the patch
              if (patchCount == iPatch) then
                 nCells = zones(iZone)%surfaceSections(surfSecCounter)%nElem
                 nConn = zones(iZone)%surfaceSections(surfSecCounter)%nConn
                 nPts = nCells*nConn
              end if
              patchCount = patchCount + 1
           end if
           surfSecCounter = surfSecCounter + 1
        end if
     end do
  end do

end subroutine getPatchSizesUnstructured

subroutine getPatchDataUnstructured(cgns_file,faceSizes,nfaces,conn,nConnPts,&
     patchPtr,patchSize,points, ndof)
  use precision
  use communication
  use gridData
  use gridInput  
  use unStructuredCGNSGrid

  implicit none
  include 'cgnslib_f.h'

  ! Input Arguments
  character*(*),intent(in) :: cgns_file
  ! Input/Ouput
  integer(kind=intType), intent(in) :: nFaces
  integer(kind=intType), intent(inout) :: faceSizes(nFaces)

  integer(kind=intType), intent(in) :: nConnPts
  integer(kind=intType), intent(inout) :: conn(nConnPts)

  integer(kind=intType), intent(in) :: patchSize
  integer(kind=intType), intent(inout) :: patchPtr(patchSize) !which nodes belong to which patch
  
  integer(kind=intType), intent(in) :: ndof
  real(kind=realType), intent(inout) :: points(ndof)

  ! Working
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
  integer(kind=intType) :: elem ,ptsIndex
  integer(kind=intType) :: elemCounter,pointCount,connCounter
  integer(kind=intType) :: patchCount
  integer(kind=intType) :: surfSecCounter


  integer(kind=intType), dimension(:,:), allocatable::localConn
  real(kind=realType), dimension(:), allocatable :: coorX, coorY, coorZ

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

     ! ! read the number of zones in the file
     ! call cg_nzones_f(cg, base, nzones, ierr)
     ! if (ierr .eq. ERROR) call cg_error_exit_f

     ! ! do a check that we have all unstructured zones
     ! do iZone = 1,nZones
     !    ! Check the zone type. This should be Unstructured. 
     !    call cg_zone_type_f(cg,base, iZone, zoneType,ierr)
     !    if (ierr .eq. ERROR) call cg_error_exit_f

     !    if(zoneType .ne. Unstructured)then
     !       print *,"Unstructured zone expected...stopping"
     !       stop
     !    end if
     ! end do

     !initialize the patch pointers
     patchPtr(1) = 0
     patchCount = 0
     elemCounter = 1
     connCounter = 1
     pointCount = 0
     ! loop over the zones and read the nodes
     zoneLoop: do iZone = 1,nZones
        surfSecCounter = 1       
        call cg_zone_read_f(cg, base, iZone, zoneName, dims, ierr)
        if (ierr .eq. CG_ERROR) call cg_error_exit_f

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

        ! Now loop over the sections in the unstructured file and
        ! figure out the boundary nodes.

        ! Determine the number of sections for this zone. Note that 
        ! surface elements can be stored in a volume zone, but they 
        ! are NOT taken into account in the number obtained from 
        ! cg_zone_read. 

        ! call cg_nsections_f(cg,base,iZone,nSections,ierr)
        ! print *, '   -> nSections',nSections    
        nSections = zones(iZone)%nSections
      
        do sec=1,nSections
           ! Determine the element type and set the pointer for the 
           ! connectivity accordingly. 
           call cg_section_read_f(cg,base,iZone,sec,secName,&
                type,eBeg,eEnd,nBdry,parentFlag,ierr)
           nElem = eEnd-eBeg+1

           call cg_npe_f(type, nconn, ierr)   


           allocate(localConn(nConn,nElem),STAT=ierr)
           localConn = 0

           ! Read the connectivity. Again, the node numbering of the 
           ! connectivities start at 1. If internally a starting
           ! index of 0 is used (typical for C-codes) 1 must be
           ! subtracted from the connectivities read. 
           
           call cg_elements_read_f(cg,base,iZone,sec,localConn,parentFlag,ierr)
           ! if (myID == 0) then
           !    print *, 'conn',shape(conn),conn(:,1)
           ! end if
           if (.not. zones(iZone)%isVolumeSection(sec)) then
              bocoType = zones(iZone)%surfaceSections(surfSecCounter)%BCType
              if (BCTypeName(bocoType) == 'BCWallViscous' .or. &
                   BCTypeName(bocoType) == 'BCWallInviscid' .or. &
                   BCTypeName(bocoType) == 'BCWall' .or. &
                   BCTypeName(bocoType) == 'BCWallViscousHeatFlux' .or. &
                   BCTypeName(bocoType) == 'BCWallViscousIsothermal') then
                 ! this is a wall surface, add the section data to the global arrays
                 do elem = 1,nElem

                    faceSizes(elemCounter)=nConn
                    elemCounter = elemCounter+1
                    do i = 1,nConn
                       conn(connCounter) = connCounter-1!pointCount!localConn(i,elem)
                       connCounter=connCounter+1
                       ptsIndex = pointCount+1
                       points(ptsIndex) = coorX(localConn(i,elem))
                       ptsIndex = pointCount+2
                       points(ptsIndex) = coorY(localConn(i,elem))
                       ptsIndex = pointCount+3
                       points(ptsIndex) = coorZ(localConn(i,elem))
                       pointCount = pointCount+3
                    end do
                 end do
                 patchPtr(patchCount+2) = elemCounter-1 !which nodes belong to which patch
                 patchCount = patchCount+1
              endif
              surfSecCounter = surfSecCounter + 1
           end if
           deallocate(localConn,STAT=ierr)
        end do
        ! And free the temporary arrays
        deallocate(coorX, coorY, coorZ)
     end do zoneLoop
  end if
end subroutine getPatchDataUnstructured
