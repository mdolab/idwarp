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
  integer(kind=intType) :: nVol, ierr, j, ii
  real(kind=realType) :: timeA, timeB, iif, err

  ! Update the nodal properties
  call computeNodalProperties(.False.)
  mytree%Ldef = Ldef0 * LdefFact
  mytree%alphaToBexp = alpha**bExp
  mytree%errTol = errTol
  nVol = size(XvPtr)/3
  call setData(mytree, Bi, Mi)
   
  timeA = mpi_wtime()
  iif = zero
  volLoop: do j=1, nVol
     r(1) = Xv0Ptr(3*j-2)
     r(2) = Xv0Ptr(3*j-1)
     r(3) = Xv0Ptr(3*j)
     num = zero
     den = zero
     ii = 0
     call evalDisp(mytree, r, num, den, ii, denomenator0(j))
     iif = iif + dble(ii)
 
     XvPtr(3*j-2) = Xv0Ptr(3*j-2) + num(1) / den
     XvPtr(3*j-1) = Xv0Ptr(3*j-1) + num(2) / den
     XvPtr(3*j  ) = Xv0Ptr(3*j  ) + num(3) / den

     ! Store these for a future sensitivity calc.
     numerator(:, j) = num
     denomenator(j) = den

  end do volLoop
  timeB = mpi_wtime()
  print *,'myid:', myid, timeB-timeA, iif/nVol, nUNIQUE

end subroutine warpMeshFast
