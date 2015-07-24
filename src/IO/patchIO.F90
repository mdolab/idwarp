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

  ! loop over the zones and read
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
  
  ! Get names one at a time since f2py can't do arrays of strings
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
  integer(kind=intType) :: patchCount
  integer(kind=intType) :: zone, boco
  integer(kind=intType) :: nZones, nBCs
  integer(kind=intType) :: bocotype

  !print *,'on Entry',iPatch!,patchName

  nZones = size(blocks) 
  patchCount = 0
  ! loop over the zones and read
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
           if (patchCount == iPatch) then
              patchName =  blocks(zone)%bocos(boco)%family
           end if
           patchCount = patchCount + 1
        end if

     end do
  end do

end subroutine getPatchNameStructured

subroutine getPatchSizeStructured(iPatch, patchSize)

  ! Get the patch sizes
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
  integer(kind=intType) :: patchCount
  integer(kind=intType) :: zone, boco
  integer(kind=intType) :: nZones, nBCs
  integer(kind=intType) :: bocotype
  integer(kind=intType) :: ptRange(3,2)

  !print *,'ipatch Size',patchIndex

  nZones = size(blocks) 
  patchCount = 0
  ! loop over the zones and read
  do zone = 1,nZones
     nBCs = size(blocks(zone)%bocos) 
     do boco=1,nBCs
        !print *,'boco num',boco,nBCs
        bocoType = blocks(zone)%bocos(boco)%bocoType
        if (BCTypeName(bocoType) == 'BCWallViscous' .or. &
             BCTypeName(bocoType) == 'BCWallInviscid' .or. &
             BCTypeName(bocoType) == 'BCWall' .or. &
             BCTypeName(bocoType) == 'BCWallViscousHeatFlux' .or. &
             BCTypeName(bocoType) == 'BCWallViscousIsothermal') then
           ! this is a wall surface, count the patch
           if (patchCount == iPatch) then
              ptRange = blocks(zone)%bocos(boco)%ptRange
              !print *,'patch size',blocks(zone)%bocos(boco)%ptRange
              if (ptRange(3,1)==ptRange(3,2)) then !Zmin/Zmax
                 patchSize(1) = blocks(zone)%bocos(boco)%ptRange(1,2)!il
                 patchSize(2) = blocks(zone)%bocos(boco)%ptRange(2,2)!jl
              else if (ptRange(1,1)==ptRange(1,2)) then ! iMin/iMax
                 patchSize(1) = blocks(zone)%bocos(boco)%ptRange(2,2)!jl
                 patchSize(2) = blocks(zone)%bocos(boco)%ptRange(3,2)!kl
              else ! jMin, jMax
                 patchSize(1) = blocks(zone)%bocos(boco)%ptRange(1,2)!il
                 patchSize(2) = blocks(zone)%bocos(boco)%ptRange(3,2)!kl
              end if
           end if
           
           patchCount = patchCount + 1
        end if

     end do
  end do
  !print *,'patchSize',patchSize
end subroutine getPatchSizeStructured

subroutine getPatchesAndConnStructured(cgns_file,points, ndof, conn,ncell)

 ! Get the patch sizes
  use precision
  use constants 
  use communication
  use structuredCGNSGrid

  implicit none
  include 'cgnslib_f.h'

 ! Input Arguments
  character*(*),intent(in) :: cgns_file
  ! Input/Ouput
  integer(kind=intType), intent(in) :: ndof
  real(kind=realType), intent(inout) :: points(ndof)

  integer(kind=intType), intent(in) :: ncell
  integer(kind=intType), intent(inout) :: conn(4*ncell)


  ! Working
  ! integer(kind=intType) :: mm, i, j, ii, ni, nj
  ! integer(kind=intType) :: zone, boco
  ! integer(kind=intType) :: nZones, nBCs
  ! integer(kind=intType) :: bocotype
  ! integer(kind=intType) :: ptRange(3,2)
  integer(kind=intType) :: i, j, k, ii, iZone
  integer(kind=intType):: ierr, base, nZones, dims(9), cg, cursize, coorsize
  integer(kind=intType):: CellDim, PhysDim
  integer(kind=intType) :: nbases, start(3), tmpSym, nSymm
  character*100 fileName
  character*32 :: zoneName, bocoName, famName, connectName, donorName, basename
  integer(kind=intType) :: npts, nbocos, bocotype, pts(3, 2), famID, boco, iB2b, transform(3)
  integer(kind=intType) :: donorRange(3, 2), nB2b
  integer(kind=intType) :: ptset_type, normalIndex(3), NormalListFlag, datatype, ndataset
  real(kind=8)   ::  data_double(6), avgNodes, symmSum(3)
  real(kind=realType), dimension(:, :, :), allocatable :: coorX, coorY, coorZ
  integer(kind=intType) :: cellCount, nodeCount, pointCount
  integer(kind=intType) :: iBeg,iEnd,jBeg,Jend,index,ni,nj
 
  ! we threw out all of the grid information when we read the CGNS file,
  ! so we need to read it a second time

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
     if (myID == 0) then
        print *, '   -> Number of Zones:', nzones
     end if

     ! the block information is already allocated, so now we just need to
     ! read the coordinates and set the surface ones into pts
     cellCount = 0
     nodeCount = 0
     pointCount = 0
     zoneLoop: do iZone=1,nZones
        !print *,'reading zones',iZone
        call cg_zone_read_f(cg, base, iZone, zoneName, dims, ierr)
        if (ierr .eq. CG_ERROR) call cg_error_exit_f
        !allocate the memory for this block and store the coords
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

        !print *,'readingbocos'
        ! get the number of boundary conditions for this zone
        call cg_nbocos_f(cg, base, iZone, nbocos, ierr)
        if (ierr .eq. ERROR) call cg_error_exit_f

        ! loop over the boundaries
        bocoLoop: do boco=1, nbocos
           !print *,'looping bocos',boco
           call cg_boco_info_f(cg, base, iZone, boco, boconame, bocotype, &
                ptset_type, npts, NormalIndex, NormalListFlag, datatype,&
                ndataset, ierr)
           if (ierr .eq. CG_ERROR) call cg_error_exit_f
           
           ! read the BC data
           call cg_boco_read_f(cg, base, iZone, boco, pts, data_double, ierr)
           if (ierr .eq. CG_ERROR) call cg_error_exit_f

           ! Determine if this is in fact a face-bc
           if ( (pts(1,1) == pts(1,2)) .or. (pts(2,1) == pts(2,2)) .or. &
                (pts(3,1) == pts(3,2))) then 
              call cg_goto_f(cg, base, ierr, "Zone_t", iZone, "ZoneBC_t", 1, "BC_t", boco, "end")
              if (ierr == 0) then ! Node exists
                 
                 ! check that it is a wall
                 if (BCTypeName(bocoType) == 'BCWallViscous' .or. &
                      BCTypeName(bocoType) == 'BCWallInviscid' .or. &
                      BCTypeName(bocoType) == 'BCWall' .or. &
                      BCTypeName(bocoType) == 'BCWallViscousHeatFlux' .or. &
                      BCTypeName(bocoType) == 'BCWallViscousIsothermal') then

                    ! set the wall nodes and generate the connectivity
                    do k=pts(3, 1), pts(3, 2)
                       do j=pts(2, 1), pts(2, 2)
                          do i=pts(1, 1), pts(1, 2)
                             index = pointCount+1
                             points(index) = coorX(i,j,k)
                             index = pointCount+2
                             points(index) = coorY(i,j,k)
                             index = pointCount+3
                             points(index) = coorZ(i,j,k)
                             pointCount = pointCount+3
                          end do
                       end do
                    end do

                    if (pts(3,1)==pts(3,2)) then !Zmin/Zmax
                       iBeg = blocks(iZone)%bocos(boco)%ptRange(1,1)!1
                       jBeg = blocks(iZone)%bocos(boco)%ptRange(2,1)!1
                       iEnd = blocks(iZone)%bocos(boco)%ptRange(1,2)!il
                       jEnd = blocks(iZone)%bocos(boco)%ptRange(2,2)!jl
                    else if (pts(1,1)==pts(1,2)) then ! iMin/iMax
                       iBeg = blocks(iZone)%bocos(boco)%ptRange(2,1)!1
                       jBeg = blocks(iZone)%bocos(boco)%ptRange(3,1)!1
                       iEnd = blocks(iZone)%bocos(boco)%ptRange(2,2)!il
                       jEnd = blocks(iZone)%bocos(boco)%ptRange(3,2)!jl
                    else ! jMin, jMax
                       iBeg = blocks(iZone)%bocos(boco)%ptRange(1,1)!1
                       jBeg = blocks(iZone)%bocos(boco)%ptRange(3,1)!1
                       iEnd = blocks(iZone)%bocos(boco)%ptRange(1,2)!il
                       jEnd = blocks(iZone)%bocos(boco)%ptRange(3,2)!jl
                    end if
           
                    ni = iEnd - iBeg + 1
                    nj = jEnd - jBeg + 1
           
                    ! Loop over generic face size...Note we are doing zero
                    ! based ordering!
                    do j=0,nj-2
                       do i=0,ni-2
                          conn(4*cellCount+1) = nodeCount + (j  )*ni + i
                          conn(4*cellCount+2) = nodeCount + (j  )*ni + i + 1
                          conn(4*cellCount+3) = nodeCount + (j+1)*ni + i + 1
                          conn(4*cellCount+4) = nodeCount + (j+1)*ni + i  
                          cellCount = cellCount + 1
                       end do
                    end do
                    nodeCount = nodeCount + ni*nj
                 end if
              end if
           end if
        end do bocoLoop

       ! free the temporary arrays
        deallocate(coorX, coorY, coorZ)
     end do zoneLoop
  end if
end subroutine getPatchesAndConnStructured
 

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

  !print *,'on Entry',iPatch!,patchName

  !nZones = size(zones) 
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
