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

  ! Determine the number of surface patches, points and cells
  nPatch = 0
  do zone = 1, size(zones)
     nPatch = nPatch + size(zones(zone)%bocos)
  end do

  ! Now we can allocate the required space
  if (.not. allocated(surfaceSizes)) then 
     allocate(surfaceSizes(2, nPatch))
     allocate(surfaceNames(nPatch))
     allocate(surfaceIsWall(nPatch))
     allocate(surfaceIsSymm(nPatch))
  end if

  surfaceIsWall = .False.
  surfaceIsSymm = .False.
  ! And make a second loop over to fill everything up!
  ! Determine the number of wall patches, points and cells
  nPatch = 0
  do zone = 1, size(zones)
     do boco=1, size(zones(zone)%bocos)
        bocoType = zones(zone)%bocos(boco)%type
        nPatch = nPatch + 1

        ! Determine if we have a wall
        if (BCTypeName(bocoType) == 'BCWallViscous' .or. &
             BCTypeName(bocoType) == 'BCWallInviscid' .or. &
             BCTypeName(bocoType) == 'BCWall' .or. &
             BCTypeName(bocoType) == 'BCWallViscousHeatFlux' .or. &
             BCTypeName(bocoType) == 'BCWallViscousIsothermal') then
           surfaceIsWall(nPatch) = .True.
        end if

        ! Determine if we have a symmetry
        if (BCTypeName(bocoType) == 'BCSymmetryPlane') then 
           surfaceIsSymm(nPatch) = .True.
        end if

        ! Set the name
        surfaceNames(nPatch) = trim(zones(zone)%bocos(boco)%family)

        ! Point range and sizes of this BC
        ptRange = zones(zone)%bocos(boco)%ptRange
        il = ptRange(1, 2) - ptRange(1, 1) + 1 ! Pts in I
        jl = ptRange(2, 2) - ptRange(2, 1) + 1 ! Pts in j
        kl = ptRange(3, 2) - ptRange(3, 1) + 1 ! Pts in K
           
        if (ptRange(1,1) ==  ptRange(1,2)) then ! iMin/iMax
           surfaceSizes(:, nPatch) = (/jl, kl/)
        else if (ptRange(2,1) == ptRange(2,2)) then ! jMin/jmax
           surfaceSizes(:, nPatch) = (/il, kl/)
        else if (ptRange(3,1) == ptRange(3,2)) then ! Kmin/Kmax
           surfaceSizes(:, nPatch) = (/il, jl/)
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
  integer(kind=intType) :: zone, boco, bocoType
  integer(kind=intType) :: nConn, nodesPerElem, nElem, cumPts
  integer(kind=intType), dimension(:), pointer :: elemPtr
  ! Determine the number of wall patches, points and cells
  nPatch = 0
  nConn = 0
  nElem =0
  do zone = 1, size(zones)
     do boco=1, size(zones(zone)%bocos)
        bocoType = zones(zone)%bocos(boco)%type
     
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
     end do
  end do

  ! Now we can allocate the required space
  if (.not. allocated(surfaceNames)) then 
     allocate(surfaceNames(nPatch))
     allocate(surfacePoints(nConn*3))
     allocate(surfaceConn(nConn))
     allocate(surfacePtr(nElem+1))
     allocate(surfacePatchPtr(nPatch+1))
     allocate(surfaceIsWall(nPatch))
     allocate(surfaceIsSymm(nPatch))
  end if

  ! Initialize the data:
  surfaceNames(:) = ""
  surfacePoints = zero
  surfaceConn = 0
  surfacePtr = 0
  surfacePatchPtr = 0
  surfaceIsWall = .False.
  surfaceIsSymm = .False.
  ! Initialize the pointers 1
  surfacePtr(1) = 1
  surfacePatchPtr(1) = 1

  ! Make a second pass back through filling things up as go.
  nPatch = 0
  nConn = 0
  cumPts = 0
  nElem = 0

  do zone = 1, size(zones)
     do boco=1, size(zones(zone)%bocos) 
        bocoType = zones(zone)%bocos(boco)%type
        nPatch = nPatch + 1

        ! Determine if we have a wall
        if (BCTypeName(bocoType) == 'BCWallViscous' .or. &
             BCTypeName(bocoType) == 'BCWallInviscid' .or. &
             BCTypeName(bocoType) == 'BCWall' .or. &
             BCTypeName(bocoType) == 'BCWallViscousHeatFlux' .or. &
             BCTypeName(bocoType) == 'BCWallViscousIsothermal') then
           surfaceIsWall(nPatch) = .True. 
        end if

        ! Determine if we have a symmetry condition
        if (BCTypeName(bocoType) == 'BCSymmetryPlane') then 
           surfaceIsSymm(nPatch) = .True.
        end if

        ! Set the name
        surfaceNames(nPatch) = trim(zones(zone)%bocos(boco)%family)

        elemPtr => zones(zone)%bocos(boco)%elemPtr
        do i=1, zones(zone)%bocos(boco)%nBCElem
           nodesPerElem = elemPtr(i+1) - elemPtr(i)
           
           do k=elemPtr(i), elemPtr(i+1)-1
              nConn = nConn + 1
              surfacePoints(3*nConn-2:3*nConn) = zones(zone)%bocos(boco)%elemNodes(:, k)
              surfaceConn(nConn) = cumPts + zones(zone)%bocos(boco)%elemConn(k) 
           end do
           nElem = nElem + 1
           surfacePtr(nElem + 1) = surfacePtr(nElem) + nodesPerElem
        end do
        cumPts = nConn
        
        ! Now many more elements did we add?
        surfacePatchPtr(nPatch + 1) = nElem + 1
     end do
  end do
 
end subroutine processUnstructuredPatches


subroutine averageNormal(pts, conn, faceSizes, nPts, nConn, nFaceSizes, avgNorm)

  use constants
  implicit none

  ! Input/Output
  real(kind=realType), intent(in) :: pts(nPts)
  integer(kind=intType), intent(in) :: conn(nConn), faceSizes(nFaceSizes)
  integer(kind=intType), intent(in) :: nPts, nConn, nFaceSizes
  real(kind=realType), intent(out) :: avgNorm(3)
  
  ! Working
  integer(kind=intType) :: i, j, nAvg, curElem, n
  real(kind=realType) :: elemPts(3,10) ! Take up to 10 sized polygon
  real(kind=realType), dimension(3) :: p1, p2, p3, a, b, c
 
  ! Determine the average normal for the patch defined by pts, conn
  ! and face sizes.

  avgNorm = zero
  nAvg = 0
  curElem = 1
  do i=1, nFaceSizes

     ! Extract out the nodes for the element
     do j=1, faceSizes(i)
        n = conn(curElem) ! This is the node number
        elemPts(:, j) = pts(3*n-2 : 3*n)
        curElem = curElem + 1
     end do
        
     ! Check triangles
     do j=0, faceSizes(i)-3

        p1 = elemPts(:, 1+j)
        p2 = elemPts(:, 2+j)
        p3 = elemPts(:, 3+j)
        
        a = p3 - p2
        b = p1 - p2

        ! Compute the cross product and normalize
        c(1) = a(2) * b(3) - a(3) * b(2)
        c(2) = a(3) * b(1) - a(1) * b(3)
        c(3) = a(1) * b(2) - a(2) * b(1)

        c = c / sqrt(c(1)**2 + c(2)**2 + c(3)**2)

        ! Now add to the 
        avgNorm = avgNorm + c
        nAvg = nAvg + 1
     end do
  end do
  avgNorm = avgNorm / nAvg

  ! Renormalize to be sure
  avgNorm = avgNorm / sqrt(avgNorm(1)**2 + avgNorm(2)**2 + avgNorm(3)**2)

end subroutine averageNormal

subroutine deallocatePatchIO
  
  use CGNSGrid

  implicit none
  
     ! Delete the surface-stuff in the cgns grid
     if (allocated(surfacePoints)) then 
        deallocate(surfacePoints)
     end if

     if (allocated(surfaceConn)) then 
        deallocate(surfaceConn)
     end if

     if (allocated(surfaceSizes)) then 
        deallocate(surfaceSizes)
     end if

     if (allocated(surfacePtr)) then 
        deallocate(surfacePtr)
     end if

     if (allocated(surfacePatchPtr)) then 
        deallocate(surfacePatchPtr)
     end if

     if (allocated(surfaceNames)) then 
        deallocate(surfaceNames)
     end if


     if (allocated(surfaceIsWall)) then 
        deallocate(surfaceIsWall)
     end if

     if (allocated(surfaceIsSymm)) then 
        deallocate(surfaceIsSymm)
     end if


end subroutine deallocatePatchIO
