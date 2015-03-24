subroutine warpDeriv(dXv_f, ndof_warp)

  use gridData
  use gridInput
  use communication
  use kd_tree
  implicit none
  
  ! Input
  integer(kind=intType), intent(in) :: ndof_warp
  real(kind=realType), intent(in) :: dXv_f(ndof_warp)
#ifndef USE_COMPLEX
  ! Working Data
  integer(kind=intType) :: i, j, kk, istart, iend, ierr, isize, nVol, nLoop
  real(kind=realType), dimension(:), allocatable :: dxsSummed
  real(kind=realType), pointer, dimension(:, :) :: Bib
  real(kind=realType), pointer, dimension(:, :, :) :: Mib
  real(kind=realType), dimension(3) :: r, rr, num, numb
  real(kind=realType) :: fact(3, 2), oden
  ! Scatter Xs into our local vector  
  call VecScatterBegin(XsToXsLocal, Xs, XsLocal, INSERT_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecScatterEnd(XsToXsLocal, Xs, XsLocal, INSERT_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Extract a pointer from Xs to use in the main routine
  call VecGetArrayF90(XsLocal, XsPtr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Extract the pointers for the volume nodes that we wish to operate on
  call VecGetArrayF90(Xv0, Xv0Ptr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecGetArrayF90(Xv, XvPtr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecGetArrayF90(dXs, dXsPtr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Allocate the extra data we need for the warping derivative:
  allocate(XvPtrb(size(XvPtr)), XsPtrb(SIZE(XsPtr)))
 
  ! dXv_f is the actual reverse seed so copy:
  XvPtrb(:) = dXv_f

  ! Zero the output surface derivative vector dXs
  call vecZeroEntries(dXs, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Now run the actual reverse mode routine
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
  Bib => mytrees(1)%tp%Bib
  Mib => mytrees(1)%tp%Mib
  call getLoopFact(nLoop, fact)
  do kk=1, nLoop
     do j=1, nvol
        r(1) = xv0ptr(3*j-2)*fact(1, kk)
        r(2) = xv0ptr(3*j-1)*fact(2, kk)
        r(3) = xv0ptr(3*j  )*fact(3, kk)

        ! Numerator derivative
        oden = one/denomenator(j)
        numb(1) = xvptrb(3*j-2)*oden*fact(1, kk)
        numb(2) = xvptrb(3*j-1)*oden*fact(2, kk)
        numb(3) = xvptrb(3*j  )*oden*fact(3, kk)
        if (evalMode == EVAL_EXACT) then 
           call evalNodeExact_b(mytrees(1)%tp, mytrees(1)%tp%root, r, numb, Bib, Mib)
        else           
           call evalNode_b(mytrees(1)%tp, mytrees(1)%tp%root, r, numb, &
                denomenator0(j), Bib, Mib)
        end if
     end do
  end do

  ! Accumulate the nodal data into Bib and Mib
  call setData_b(mytrees(1)%tp, mytrees(1)%tp%root)

  ! And now accumulate into XsPtrb
  call computeNodalProperties_b(mytrees(1)%tp, .False.)
  
  ! And remove our derviative values
  call deallocDerivValues(mytrees(1)%tp)

  call VecGetSize(dXs, isize, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecGetOwnershipRange(Xs, istart, iend, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  allocate(dxssummed(isize))
  dxssummed = zero

  call MPI_Allreduce(xsptrb, dxssummed, isize, MPI_REAL8, MPI_SUM, warp_comm_world, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  
  ! Copy just the required part into dXsPtr
  dXsPtr = dxssummed(iStart+1:iEnd)
  deallocate(dxssummed)
  
  ! Deallocate the extra space
  deallocate(XvPtrb, XsPtrb)

  ! Restore all the arrays
  call VecRestoreArrayF90(XsLocal, XsPtr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecRestoreArrayF90(Xv0, Xv0Ptr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecRestoreArrayF90(Xv, XvPtr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecRestoreArrayF90(dXs, dXsPtr, ierr)
  call EChk(ierr,__FILE__,__LINE__)
#endif  
end subroutine warpDeriv
