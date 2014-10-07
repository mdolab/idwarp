subroutine warpMeshFast

  ! This is gateway function to the NlogN (fast) mesh warping
  ! algorithm. The input to this function is XsPtr a (flatted) 1D
  ! array of the surface coordinates and the output is XvPtr a
  ! (flatted) 1D arary of the volume coordinates.

  use gridData
  use gridInput
  use communication
  use kd_tree
  implicit none

  ! Working parameters
  real(kind=realType) :: den
  real(kind=realType), dimension(3) :: r, num
  integer(kind=intType) :: nVol, ierr, j, kk, ii, nLoop
  real(kind=realType) :: fact(3, 2)
  denomenator = zero
  numerator = zero
  nVol = size(XvPtr)/3
  call computeNodalProperties(mytrees(1)%tp, .False.)
  mytrees(1)%tp%Ldef = mytrees(1)%tp%Ldef0 * LdefFact
  mytrees(1)%tp%alphaToBexp = alpha**bExp
  mytrees(1)%tp%errTol = errTol
  call setData(mytrees(1)%tp, mytrees(1)%tp%root)
  call getLoopFact(nLoop, fact)

  do kk=1, nLoop
     volLoop: do j=1, nVol
        r(1) = Xv0Ptr(3*j-2)*fact(1, kk)
        r(2) = Xv0Ptr(3*j-1)*fact(2, kk)
        r(3) = Xv0Ptr(3*j)*fact(3, kk)

        num = zero
        den = zero 
        call evalDisp(mytrees(1)%tp, mytrees(1)%tp%root, r, num, den, denomenator0(j))

        numerator(1, j) = numerator(1, j) + num(1)*fact(1, kk)
        numerator(2, j) = numerator(2, j) + num(2)*fact(2, kk)
        numerator(3, j) = numerator(3, j) + num(3)*fact(3, kk)
        denomenator(j) = denomenator(j) + den
     end do volLoop
  end do

  updateLoop: do j=1, nVol
     XvPtr(3*j-2) = Xv0Ptr(3*j-2) + numerator(1, j) / denomenator(j)
     XvPtr(3*j-1) = Xv0Ptr(3*j-1) + numerator(2, j) / denomenator(j)
     XvPtr(3*j  ) = Xv0Ptr(3*j  ) + numerator(3, j) / denomenator(j)
  end do updateLoop
end subroutine warpMeshFast
