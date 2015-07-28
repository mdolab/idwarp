! **************************************************
! Routines for dealing with the structured CGNS mesh
! **************************************************

subroutine processStructuredPatches

  ! This subroutine does the required preprocessing that converts all
  ! the bounary condition information into a connected surface mesh
  ! suitable for use in the unstructured mesh warping algorithm.

  use CGNSGrid

  implicit none
  include 'cgnslib_f.h'

  ! Working variables
  integer(kind=intType) :: nPatch
  integer(kind=intType) :: zone, boco, bocoType, ptRange(3, 2)
  integer(kind=intType) :: il, jl, kl

  ! Determine the number of wall patches, points and cells
  nPatch = 0
  do zone = 1, size(zones)
     do boco=1, size(zones(zone)%bocos) 
        bocoType = zones(zone)%bocos(boco)%type
        if (BCTypeName(bocoType) == 'BCWallViscous' .or. &
             BCTypeName(bocoType) == 'BCWallInviscid' .or. &
             BCTypeName(bocoType) == 'BCWall' .or. &
             BCTypeName(bocoType) == 'BCWallViscousHeatFlux' .or. &
             BCTypeName(bocoType) == 'BCWallViscousIsothermal') then

           ! this is a wall surface, count the patch
           nPatch = nPatch + 1
        end if
     end do
  end do

  ! Now we can allocate the required space
  allocate(wallSizes(2, nPatch))
  allocate(wallNames(nPatch))

  ! And make a second loop over to fill everything up!
  ! Determine the number of wall patches, points and cells
  nPatch = 0
  do zone = 1, size(zones)
     do boco=1, size(zones(zone)%bocos)
        bocoType = zones(zone)%bocos(boco)%type
        if (BCTypeName(bocoType) == 'BCWallViscous' .or. &
             BCTypeName(bocoType) == 'BCWallInviscid' .or. &
             BCTypeName(bocoType) == 'BCWall' .or. &
             BCTypeName(bocoType) == 'BCWallViscousHeatFlux' .or. &
             BCTypeName(bocoType) == 'BCWallViscousIsothermal') then

           ! this is a wall surface, count the patch
           nPatch = nPatch + 1

           ! Set the name
           wallNames(nPatch) = trim(zones(zone)%bocos(boco)%family)

           ! And the sizes:

           ! Point range and sizes of this BC
           ptRange = zones(zone)%bocos(boco)%ptRange
           il = ptRange(1, 2) - ptRange(1, 1) + 1 ! Pts in I
           jl = ptRange(2, 2) - ptRange(2, 1) + 1 ! Pts in j
           kl = ptRange(3, 2) - ptRange(3, 1) + 1 ! Pts in K
           
           if (ptRange(1,1) ==  ptRange(1,2)) then ! iMin/iMax
              wallSizes(:, nPatch) = (/jl, kl/)
           else if (ptRange(2,1) == ptRange(2,2)) then ! jMin/jmax
              wallSizes(:, nPatch) = (/il, kl/)
           else if (ptRange(3,1) == ptRange(3,2)) then ! Kmin/Kmax
              wallSizes(:, nPatch) = (/il, jl/)
           end if
        end if
     end do
  end do
end subroutine processStructuredPatches

subroutine processUnstructuredPatches
  
  ! This subroutine does the required preprocessing that converts all
  ! the boundary condition information into a connected surface mesh
  ! suitable for use in the unstructured mesh warping algorithm. We
  ! essentially need to create cgnsGrid.wallConn, cgnsGrid.wallPts,
  ! cgnsGrid.patchName and wallPtr. Unlike the structured code above,
  ! wallPts and wallConn are not yet allocated and set. 

  use CGNSGrid

  implicit none
  include 'cgnslib_f.h'

  ! Working variables
  integer(kind=intType) :: nPatch, i, k
  integer(kind=intType) :: zone, boco, bocoType, ptRange(3, 2)
  integer(kind=intType) :: nConn, nodesPerElem, nElem, cumPts
  integer(kind=intType), dimension(:), pointer :: elemPtr
  ! Determine the number of wall patches, points and cells
  nPatch = 0
  nConn = 0
  nElem =0
  do zone = 1, size(zones)
     do boco=1, size(zones(zone)%bocos)
        bocoType = zones(zone)%bocos(boco)%type
        if (BCTypeName(bocoType) == 'BCWallViscous' .or. &
             BCTypeName(bocoType) == 'BCWallInviscid' .or. &
             BCTypeName(bocoType) == 'BCWall' .or. &
             BCTypeName(bocoType) == 'BCWallViscousHeatFlux' .or. &
             BCTypeName(bocoType) == 'BCWallViscousIsothermal') then

           ! this is a wall surface, so count the patch
           nPatch = nPatch + 1
           
           ! Now determine how many points we will have and the size
           ! we need need for the connectivity

           ! We have to be a little careful here since nBCElem is the
           ! *actual* number of BCElements and may be less than size(BCElements).
           elemPtr => zones(zone)%bocos(boco)%elemPtr
           do i=1, zones(zone)%bocos(boco)%nBCElem
              nodesPerElem = elemPtr(i+1) - elemPtr(i)
              nConn = nConn + nodesPerElem
           end do
           nElem = nElem + zones(zone)%bocos(boco)%nBCElem
        end if
     end do
  end do

  ! Now we can allocate the required space
  allocate(wallNames(nPatch))
  allocate(wallPoints(nConn*3))
  allocate(wallConn(nConn))
  allocate(wallPtr(nElem+1))
  allocate(wallPatchPtr(nPatch+1))
  
  ! Initialize the data:
  wallNames(:) = ""
  wallPoints = zero
  wallConn = 0
  wallPtr = 0
  wallPatchPtr = 0

  ! Initialize the pointers 1
  wallPtr(1) = 1
  wallPatchPtr(1) = 1

  ! Make a second pass back through filling things up as go.
  nPatch = 0
  nConn = 0
  cumPts = 0
  nElem = 0
  do zone = 1, size(zones)
     do boco=1, size(zones(zone)%bocos) 
        bocoType = zones(zone)%bocos(boco)%type
        if (BCTypeName(bocoType) == 'BCWallViscous' .or. &
             BCTypeName(bocoType) == 'BCWallInviscid' .or. &
             BCTypeName(bocoType) == 'BCWall' .or. &
             BCTypeName(bocoType) == 'BCWallViscousHeatFlux' .or. &
             BCTypeName(bocoType) == 'BCWallViscousIsothermal') then

           ! this is a wall surface, so count the patch
           nPatch = nPatch + 1
           
           ! Set the name
           wallNames(nPatch) = trim(zones(zone)%bocos(boco)%family)

           elemPtr => zones(zone)%bocos(boco)%elemPtr
           do i=1, zones(zone)%bocos(boco)%nBCElem
              nodesPerElem = elemPtr(i+1) - elemPtr(i)
              
              do k=elemPtr(i), elemPtr(i+1)-1
                 nConn = nConn + 1
                 wallPoints(3*nConn-2:3*nConn) = zones(zone)%bocos(boco)%elemNodes(:, k)
                 wallConn(nConn) = cumPts + zones(zone)%bocos(boco)%elemConn(k) 
              end do
              nElem = nElem + 1
              wallPtr(nElem + 1) = wallPtr(nElem) + nodesPerElem
           end do
           cumPts = nConn

           ! Now many more elements did we add?
           wallPatchPtr(nPatch + 1) = cumPts+1
        end if
     end do
  end do
 
end subroutine processUnstructuredPatches
