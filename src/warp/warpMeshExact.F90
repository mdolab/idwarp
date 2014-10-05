subroutine warpMeshExact

  ! This is gateway function to the N^2 (exact) mesh warping
  ! algorithm. The input to this function is XsPtr a (flatted) 1D
  ! array of the surface coordinates and the ouput XvPtr a (flatted)
  ! 1D arary of the volume coordinates. This routine is AD-able by
  ! Tapendade in both forward and reverse modes.

  use gridData
  use gridInput
  use communication
  implicit none

  ! Working parameters
  real(kind=realType) :: dist, LdefoDist, Wi, den
  real(kind=realType), dimension(3) :: r, rr,  num, Si
  integer(kind=intType) :: nVol, i, j
  real(kind=realType) :: Ldef, alphaToBexp, oDen
  real(kind=realType) :: LdefoDist3
  real(kind=realType) :: timeA, timeB

  ! Update the nodal properties
  call computeNodalProperties(.False.)

  ! Compute the actual Ldef
  Ldef = Ldef0 * LdefFact
  
  ! Pre-compute this
  alphaToBexp = alpha**bExp
  
  ! Allocate as needed
  nVol = size(XvPtr)/3

  numerator= zero
  denomenator = zero
  ! Loop over the volume nodes
#ifndef USE_TAPENADE
   timeA = mpi_wtime()
#endif

  !$AD II-LOOP
  volLoop: do j=1, nVol
     r(1) = Xv0Ptr(3*j-2)
     r(2) = Xv0Ptr(3*j-1)
     r(3) = Xv0Ptr(3*j)
     num = zero
     den = zero
     if (useRotations) then 
        do i=1, nUnique
           rr = r - Xu0(:, i)
           dist = sqrt(rr(1)**2 + rr(2)**2 + rr(3)**2+1e-16)
           LdefoDist = Ldef/dist
           Ldefodist3 = LdefoDist**3
           Wi = Ai(i)*(Ldefodist3 + alphaToBexp*Ldefodist3*LdefoDist*LdefoDist)
           Si = Mi(:, 1, i)*r(1) + Mi(:, 2, i)*r(2) + Mi(:, 3, i)*r(3) + bi(:, i) - r(:)
           num = num + Wi*Si
           den = den + Wi
        end do
     else
        do i=1, nUnique
           rr = r - Xu0(:, i)
           dist = sqrt(rr(1)**2 + rr(2)**2 + rr(3)**2 + 1e-15)
           LdefoDist = Ldef/dist
           Ldefodist3 = LdefoDist**3
           Wi = Ai(i)*(Ldefodist3 + alphaToBexp*Ldefodist3*LdefoDist*LdefoDist)
           num = num + Wi*Bi(:, i)
           den = den + Wi
        end do
     end if

     denomenator(j) = den
     numerator(:, j) = num

  end do volLoop

  ! We can do the actual update
  !$AD II-LOOP
  do j=1, nVol
     oDen = one/denomenator(j)
     XvPtr(3*j-2) = Xv0Ptr(3*j-2) + numerator(1, j) * oDen
     XvPtr(3*j-1) = Xv0Ptr(3*j-1) + numerator(2, j) * oDen
     XvPtr(3*j  ) = Xv0Ptr(3*j  ) + numerator(3, j) * oDen
  end do

end subroutine warpMeshExact


subroutine getCaseFromExponents(aExp, bExp, iCase)
  use constants
  implicit none
  real(kind=realType), intent(in) :: aExp, bExp
  integer(kind=intType), intent(out) :: iCase
  real(kind=realType) :: expTol
  expTol = 0.01_realType
  iCase = 0

  if      (abs(aExp - two) < expTol .and. abs(bExp - two  ) < expTol) then 
     iCase = 1
  else if (abs(aExp - two) < expTol .and. abs(bExp - three) < expTol) then 
     iCase = 2
  else if (abs(aExp - two) < expTol .and. abs(bExp - four ) < expTol) then 
     iCase = 3
  else if (abs(aExp - two) < expTol .and. abs(bExp - five ) < expTol) then 
     iCase = 4

  else if (abs(aExp - three) < expTol .and. abs(bExp - two  ) < expTol) then 
     iCase = 5
  else if (abs(aExp - three) < expTol .and. abs(bExp - three) < expTol) then 
     iCase = 6
  else if (abs(aExp - three) < expTol .and. abs(bExp - four ) < expTol) then 
     iCase = 7
  else if (abs(aExp - three) < expTol .and. abs(bExp - five ) < expTol) then 
     iCase = 8

  else if (abs(aExp - four) < expTol .and. abs(bExp - two  ) < expTol) then 
     iCase = 9
  else if (abs(aExp - four) < expTol .and. abs(bExp - three) < expTol) then 
     iCase = 10
  else if (abs(aExp - four) < expTol .and. abs(bExp - four ) < expTol) then 
     iCase = 11
  else if (abs(aExp - four) < expTol .and. abs(bExp - five ) < expTol) then 
     iCase = 12
  end if

end subroutine getCaseFromExponents
