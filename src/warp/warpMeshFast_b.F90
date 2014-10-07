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
  real(kind=realType), dimension(3) :: r, rr, num, numb
  integer(kind=intType) :: nVol, ierr, j, i, kk, nSurf, nLoop
  real(kind=realType), pointer, dimension(:, :) :: Xu0, Bib
  real(kind=realType) :: fact(3, 2)
  nVol = size(XvPtr)/3
  XsPtrb = zero
  call computeNodalProperties(mytrees(1)%tp, .False.)
  mytrees(1)%tp%Ldef = mytrees(1)%tp%Ldef0 * LdefFact
  mytrees(1)%tp%alphaToBexp = alpha**bExp
  mytrees(1)%tp%errTol = errTol

  call setData(mytrees(1)%tp, mytrees(1)%tp%root)
  call allocDerivValues(mytrees(1)%tp)
  call zeroDeriv(mytrees(1)%tp%root)

  ! Extract pointers for easier reading:
  Xu0 => mytrees(1)%tp%Xu0
  Bib => mytrees(1)%tp%Bib
  nSurf = mytrees(1)%tp%n

  call getLoopFact(nLoop, fact)
  Bib = zero
  do kk=1, nLoop
     do j=1, nvol
        r(1) = xv0ptr(3*j-2)*fact(1, kk)
        r(2) = xv0ptr(3*j-1)*fact(2, kk)
        r(3) = xv0ptr(3*j  )*fact(3, kk)
        den = denomenator(j)

        ! Numerator derivative
        numb(1) = xvptrb(3*j-2)/den
        numb(2) = xvptrb(3*j-1)/den
        numb(3) = xvptrb(3*j  )/den

        call evaLDisp_b(mytrees(1)%tp, mytrees(1)%tp%root, r, numb, denomenator0(j), &
             Bib(:, :), fact(:, kk))
     end do
  end do

  ! Accumulate the nodal data into Bib
  call setData_b(mytrees(1)%tp, mytrees(1)%tp%root)

  ! And now accumulate into XsPtrb
  call computeNodalProperties_b(mytrees(1)%tp, .False.)
  
  ! And remove our derviative values
  call deallocDerivValues(mytrees(1)%tp)
END SUBROUTINE warpMeshFast_b
