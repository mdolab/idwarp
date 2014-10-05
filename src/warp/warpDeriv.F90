subroutine warpDeriv(dXv_f, ndof_warp)

  use gridData
  use gridInput
  implicit none

  ! Input
  integer(kind=intType), intent(in) :: ndof_warp
  real(kind=realType), intent(in) :: dXv_f(ndof_warp)

  ! Working Data
  integer(kind=intType) :: i, iend, ierr

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

  ! Allocate the extra data we need for the warping derivative:
  allocate(XvPtrb(size(XvPtr)), XsPtrb(SIZE(XsPtr)))
  allocate(xub(3, nUnique))
  allocate(Mib(3, 3, nUnique), Bib(3, nUnique))
  allocate(normalsb(3, nUnique), normals0b(3, nUnique))
 
  ! dXv_f is the actual reverse seed so copy:
  XvPtrb(:) = dXv_f

  ! Zero the output surface derivative vector dXs
  call vecZeroEntries(dXs, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Now run the actual reverse mode routine
  if (evalMode == EVAL_EXACT) then 
     call warpMeshExact_b()
  else
     call warpMeshFast_b()
  end if

  ! We add ALL values from Xsptb on all procs. 
  call VecGetSize(dXs, iend, ierr)
  do i=0, iend-1
     call VecSetValues(dXs, 1, (/i/), (/XsPtrb(i+1)/), ADD_VALUES, ierr)
     call EChk(ierr, __FILE__, __LINE__)
  end do
  
  ! Must assemble 
  call VecAssemblyBegin(dXs, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecAssemblyEnd(dXs, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Deallocate the extra space
  deallocate(XvPtrb, XsPtrb)
  deallocate(mib, Bib, normalsb, normals0b, xub)

  ! Restore all the arrays
  call VecRestoreArrayF90(XsLocal, XsPtr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecRestoreArrayF90(Xv0, Xv0Ptr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecRestoreArrayF90(Xv, XvPtr, ierr)
  call EChk(ierr,__FILE__,__LINE__)
end subroutine warpDeriv
