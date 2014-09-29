subroutine computeNodalProperties(initialPoint)

  use gridData
  use gridInput
  use communication
  implicit none

  ! Subroutine Variable
  logical :: initialPoint

  ! Working variables
  real(kind=realType), allocatable, dimension(:, :) :: points
  integer(kind=intType) :: i, j, ierr, nPts, nElem, ind
  integer(kind=intType) :: nPtsMax = 10
  real(kind=realType), dimension(:, :), allocatable :: faceNormals
  real(kind=realType), dimension(:), allocatable :: faceAreas
  real(kind=realType), pointer, dimension(:) :: xx
  real(kind=realType) :: sumArea, sumNormal(3), Si(3), ds(3), sMean(3), da, eta, r(3), dx(3)

  ! Scatter Xs into our local vector  
  call VecScatterBegin(XsToXsLocal, Xs, XsLocal, INSERT_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecScatterEnd(XsToXsLocal, Xs, XsLocal, INSERT_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Extract a pointer from XsLocal and full up Xu
  call VecGetArrayF90(XsLocal, xx, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! This performs the copy and mirroring as required
  do i=1, nUnique
     j = XuInd(i)
     Xu(:, i) = xx(3*j-2:3*j)*XuFact(:, i)
  end do
  call VecRestoreArrayF90(XsLocal, xx, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! First we loop over each of the surface elements and compute the
  ! required properties:  the area and the normal vector
  allocate(points(3, nPtsMax))
  allocate(faceAreas(nFace), faceNormals(3, nFace))
  faceAreas = zero
  faceNormals = zero

  do i=1, nFace
     nPts = facePtr(i) - facePtr(i-1)

     ! Check if we need to realloc
     if (nPts > nPtsMax) then
        deallocate(points)
        nPtsMax = nPtsMax * 2
        allocate(points(3, nPtsMax))
     end if
     ! Copy out the points
     do j=1, nPts
        points(:, j) = Xu(:, faceConn(facePtr(i-1) + j))
     end do

     if (i==1) then 
        call getElementProps(points, nPts, faceAreas(i), faceNormals(:, i))
     else
        call getElementProps(points, nPts, faceAreas(i), faceNormals(:, i))
     end if
  end do

  ! Now we loop over nodes and fill up Ai, Bi, and Mi
  normals = zero
  mi = zero
  bi = zero

  do i=1, nUnique
     sumArea = zero
     sumNormal = zero
     nElem = nodeToElem(1, i)

     do j=1, nElem
        ind = nodeToElem(1+j, i)
        ! For face 'ind' how many nodes are on the element?
        dA = faceAreas(ind) / (facePtr(ind) - facePtr(ind-1))
        sumArea = sumArea + dA
        sumNormal = sumNormal + dA*facenormals(:, ind)
     end do
     Ai(i) = sumArea
     normals(:, i) = sumNormal / sumArea
     if (initialPoint) then
        normals0(:, i) = normals(:, i)
     else
   
        ! Now get the rotation Matrix
        call getRotationMatrix3d(normals0(:, i), normals(:, i), Mi(:, :, i))
           
        ! Actual calc....expand matmul for efficiency:
        !Bi(:, i) = Xu(:, i) - matmul(Mi(:, :, i), Xu0(:, i))
        Bi(:, i) = Xu(:, i) - (Mi(:, 1, i)*Xu0(1, i) + Mi(:, 2, i)*Xu0(2, i) + Mi(:, 3, i)*Xu0(3, i))
     end if
  end do

  deallocate(faceAreas, faceNormals)

  ! -------------------------------------------------------
  !        Automatic Alpha calc...not used becuase of 
  !        derivative dependancy issues
  ! -------------------------------------------------------
  ! smean = zero
  ! do i=1,nUnique
  !    r = Xu(:, i)
  !    Si(1) = Mi(1, 1, i)*r(1) + Mi(1, 2, i)*r(2) + Mi(1, 3, i)*r(3) + bi(1, i) - r(1)
  !    Si(2) = Mi(2, 1, i)*r(1) + Mi(2, 2, i)*r(2) + Mi(2, 3, i)*r(3) + bi(2, i) - r(2)
  !    Si(3) = Mi(3, 1, i)*r(1) + Mi(3, 2, i)*r(2) + Mi(3, 3, i)*r(3) + bi(3, i) - r(3)
  !    smean = smean + Ai(i) * Si
  ! end do
  ! smean = smean / sum(Ai)

  ! eta = 5.0
  ! alpha = zero
  ! do i=1, nUnique
  !    r = Xu(:, i)
  !    Si(1) = Mi(1, 1, i)*r(1) + Mi(1, 2, i)*r(2) + Mi(1, 3, i)*r(3) + bi(1, i) - r(1)
  !    Si(2) = Mi(2, 1, i)*r(1) + Mi(2, 2, i)*r(2) + Mi(2, 3, i)*r(3) + bi(2, i) - r(2)
  !    Si(3) = Mi(3, 1, i)*r(1) + Mi(3, 2, i)*r(2) + Mi(3, 3, i)*r(3) + bi(3, i) - r(3)

  !    ds = Si - sMean
  !    alpha = max(alpha, sqrt(ds(1)**2 + ds(2)**2 + ds(3)**2))
  ! end do
  ! !alpha = alpha * eta / Ldef



end subroutine computeNodalProperties

