! This is a hand-coded derivative routine for warpMeshFast. It draws
! upon the tapenade reverse AD code for the exact code. However, it is
! completely written by hand and cannot be automatically generated. 

subroutine warpMeshFast_b()
  use gridInput
  USE gridData
  use kd_tree
  IMPLICIT NONE

  ! Working parameters
  real(kind=realType) :: den
  real(kind=realType), dimension(3) :: r, num, numb
  integer(kind=intType) :: nVol, ierr, j, i

  ! Update the nodal properties
  call computeNodalProperties(.False.)
  mytree%Ldef = Ldef0 * LdefFact
  mytree%alphaToBexp = alpha**bExp
  mytree%errTol = errTol
  nVol = size(XvPtr)/3
  call setData(mytree, Bi, Mi)
  Bib = zero
  mytree%Bib => Bib
  call zeroDeriv(mytree%root)
  nvol = SIZE(xvptr)/3

  DO j=1,nvol
     r(1) = xv0ptr(3*j-2)
     r(2) = xv0ptr(3*j-1)
     r(3) = xv0ptr(3*j)
     den = denomenator(j)

     ! Numerator derivative
     numb(1) = xvptrb(3*j-2)/den
     numb(2) = xvptrb(3*j-1)/den
     numb(3) = xvptrb(3*j  )/den

     call evaLDisp_b(mytree, r, numb, denomenator0(j), Bib)
  end DO
  ! Accumulate into Bib
  call setData_b(myTree)

  CALL COMPUTENODALPROPERTIES_B(.false.)
END SUBROUTINE warpMeshFast_b
