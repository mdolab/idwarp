subroutine computeNodalProperties(initialPoint)

  use gridData
  use gridInput
  use communication
  implicit none

  ! Subroutine Variable
  logical :: initialPoint

  ! Working variables
  real(kind=realType), dimension(3, 20) :: points
  integer(kind=intType) :: i, j, jj, kk, nPts, nElem, ind
  integer(kind=intType) :: nPtsMax = 10
  real(kind=realType) :: faceArea, faceNormal(3)
  real(kind=realType) :: sumArea, sumNormal(3), Si(3), ds(3), sMean(3), da, eta, r(3), dx(3)

  ! This performs the copy and mirroring as required
  do i=1, nUnique
     j = XuInd(i)
     Xu(:, i) = XsPtr(3*j-2:3*j)*XuFact(:, i)
  end do

  ! Now we loop over nodes and fill up Ai, Bi, and Mi
  normals = zero
  mi = zero
  bi = zero

  do i=1, nUnique
     sumArea = zero
     sumNormal = zero
     nElem = nodeToElem(1, i)

     do jj=1, nElem
        ind = nodeToElem(1+jj, i)

        ! Extract points for this face
        nPts = facePtr(ind) - facePtr(ind-1)
        do kk=1, nPts
           points(:, kk) = Xu(:, faceConn(facePtr(ind-1) + kk))
        end do

        call getElementProps(points, nPts, faceArea, faceNormal)

        ! For face 'ind' how many nodes are on the element?
        dA = faceArea / (facePtr(ind) - facePtr(ind-1))
        sumArea = sumArea + dA
        sumNormal = sumNormal + dA*facenormal
     end do
     normals(:, i) = sumNormal / sumArea

     if (initialPoint) then
        normals0(:, i) = normals(:, i)
        Ai(i) = sumArea
     else
        ! Now get the rotation Matrix
        call getRotationMatrix3d(normals0(:, i), normals(:, i), Mi(:, :, i))
        if (useRotations) then
           ! Actual calc....expand matmul for efficiency:
           Bi(:, i) = Xu(:, i) - (Mi(:, 1, i)*Xu0(1, i) + Mi(:, 2, i)*Xu0(2, i) + Mi(:, 3, i)*Xu0(3, i))
        else
           Bi(:, i) = Xu(:, i) - Xu0(:, i)
        end if
     end if
  end do

end subroutine computeNodalProperties

