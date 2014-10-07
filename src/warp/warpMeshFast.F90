subroutine warpMeshFast

  ! This is gateway function to the NlogN (fast) mesh warping
  ! algorithm. The input to this function is XsPtr a (flatted) 1D
  ! array of the surface coordinates and the output is XvPtr a
  ! (flatted) 1D arary of the volume coordinates. This routine is
  ! AD-able by Tapendade in both forward and reverse modes.

  use gridData
  use gridInput
  use communication
  use kd_tree
  implicit none

  ! Working parameters
  real(kind=realType) :: den
  real(kind=realType), dimension(3) :: r, num
  integer(kind=intType) :: nVol, ierr, j, kk


  denomenator = zero
  numerator = zero
  nVol = size(XvPtr)/3

  do kk=1,size(mytrees)
     ! Update the nodal properties
     call computeNodalProperties(mytrees(kk)%tp, .False.)
     mytrees(kk)%tp%Ldef = mytrees(kk)%tp%Ldef0 * LdefFact
     mytrees(kk)%tp%alphaToBexp = alpha**bExp
     mytrees(kk)%tp%errTol = errTol

     call setData(mytrees(kk)%tp, mytrees(kk)%tp%root)

     volLoop: do j=1, nVol
        r(1) = Xv0Ptr(3*j-2)
        r(2) = Xv0Ptr(3*j-1)
        r(3) = Xv0Ptr(3*j)
        num = zero
        den = zero 
        call evalDisp(mytrees(kk)%tp, mytrees(kk)%tp%root, r, num, den, denomenator0(j))

        numerator(1, j) = numerator(1, j) + num(1)
        numerator(2, j) = numerator(2, j) + num(2)
        numerator(3, j) = numerator(3, j) + num(3)
        denomenator(j) = denomenator(j) + den
     end do volLoop
  end do

  updateLoop: do j=1, nVol
     XvPtr(3*j-2) = Xv0Ptr(3*j-2) + numerator(1, j) / denomenator(j)
     XvPtr(3*j-1) = Xv0Ptr(3*j-1) + numerator(2, j) / denomenator(j)
     XvPtr(3*j  ) = Xv0Ptr(3*j  ) + numerator(3, j) / denomenator(j)
  end do updateLoop
  
end subroutine warpMeshFast
