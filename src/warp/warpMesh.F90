subroutine warpMesh()

  ! This is the master mesh warping computation. It is necessary that
  ! comptueNodalProperties is called before this routine is called.

  use gridData
  use gridInput
  use communication
  implicit none

  ! Working parameters
  real(kind=realType), pointer, dimension(:) :: xxv0, xxv
  real(kind=realType) :: dist, LdefoDist, Wi, Si(3), denomenator
  real(kind=realType), dimension(3) :: r, rr, numerator, dx
  integer(kind=intType) :: nVol, ierr, i, j, combo
  real(kind=realType) :: timeA, timeB, Ldef, alphaToBexp, Ldefodist2, Ldefodist3, Ldefodist4, Ldefodist5
  real(kind=realType) :: intPow
  !real(kind=realType) :: aaexp, bbexp
  ! Extract the pointer for the volume nodes that we wish to operate on
  call VecGetArrayF90(Xv0, xxv0, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecGetArrayF90(Xv, xxv, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecGetLocalSize(Xv0, nVol, ierr)
  call EChk(ierr,__FILE__,__LINE__)
  nVol = nVol / 3 

  ! Compute the actual Ldef
  Ldef = Ldef0 * LdefFact

  timeA = mpi_wtime()
  alphaToBexp = alpha**bExp

  if (aexp == 4 .and. bexp == 5) then
      combo = 2
   else if (aexp ==3 .and. bexp == 5) then
      combo = 1
   else
      print *, 'bad combo'
      stop
   end if
   print *, 'combo', combo
   !aaexp = 3.2345345
   !bbexp = 4.5335464256
  ! Loop over the volume nodes
  volLoop: do j=1, nVol
     if (mod(j, 10000) == 0) then 
        print *, myid, j
     end if
     r = xxv0(3*j-2:3*j)
     numerator= zero
     denomenator = zero
     surfLoop: do i=1, nUnique
        rr = r - Xu0(:, i)
        dist = sqrt(rr(1)**2 + rr(2)**2 + rr(3)**2+1e-16)

        LdefoDist = Ldef/dist


        select case(combo)
        case (1)
           Ldefodist3 = LdefoDist**3
           Ldefodist5 = LdefoDist3*Ldefodist*Ldefodist
           Wi = Ldefodist3  + alphaToBexp*Ldefodist5
        case (2)
           Ldefodist4 = LdefoDist**4
           Ldefodist5 = LdefoDist4*Ldefodist
           Wi = Ldefodist4  + alphaToBexp*Ldefodist5
        end select
        !  Wi = Ldefodist**aaExp + (alpha*Ldefodist)**bexp

        Wi = Wi * Ai(i)

        ! Si(1) = Mi(1, 1, i)*r(1) + Mi(1, 2, i)*r(2) + Mi(1, 3, i)*r(3) + bi(1, i) - r(1)
        ! Si(2) = Mi(2, 1, i)*r(1) + Mi(2, 2, i)*r(2) + Mi(2, 3, i)*r(3) + bi(2, i) - r(2)
        ! Si(3) = Mi(3, 1, i)*r(1) + Mi(3, 2, i)*r(2) + Mi(3, 3, i)*r(3) + bi(3, i) - r(3)

        Si = Mi(:, 1, i)*r(1) + Mi(:, 2, i)*r(2) + Mi(:, 3, i)*r(3) + bi(:, i) - r(:)

        numerator = numerator + Wi*Si
        denomenator = denomenator + Wi
     end do surfLoop
     dx = numerator/denomenator

     ! Actual update
     xxv(3*j-2:3*j) = xxv0(3*j-2:3*j) + dx
  end do volLoop

  call VecRestoreArrayF90(Xv0, xxv0, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecRestoreArrayF90(Xv, xxv, ierr)
  call EChk(ierr,__FILE__,__LINE__)
  timeB = mpi_wtime()

  print *,' Mesh warp Time:', timeB-timeA
end subroutine warpMesh

