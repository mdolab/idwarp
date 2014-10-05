subroutine warpDeriv(dXv_f, ndof_warp)

  use gridData
  use communication
  use diffSizes
  implicit none

  ! Input
  integer(kind=intType), intent(in) :: ndof_warp
  real(kind=realType), intent(in) :: dXv_f(ndof_warp)

  ! Working Data
  integer(kind=intType) :: i, istart, iend, ierr
  real(kind=realType) :: timeA, timeB
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
  allocate(Mib(3, 3, nUnique), Bib(3, nUnique), Aib(nUnique))
  allocate(normalsb(3, nUnique), normals0b(3, nUnique))
  allocate(numeratorb(3, size(XvPtr)/3), denomenatorb(size(XvPtr)))
  isize1ofdrfnormals = 3
  isize2ofdrfnormals = nUnique

  isize1ofdrfnormals0 = 3
  isize2ofdrfnormals0 = nUnique

  isize1ofdrfxu = 3
  isize2ofdrfxu = nUnique

  ! dXv_f is the actual reverse seed so copy:
  XvPtrb(:) = dXv_f

  ! Zero the output surface derivative vector dXs
  call vecZeroEntries(dXs, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  print *,'calling reverse mode code...'
  ! Now run the actual reverse mode routine
  timeA = mpi_wtime()
  call warpMeshExact_b()
  timeB = mpi_wtime()
  print *,'done calling reverse mode code', timeB - timeA
  ! And we now need take local values of XsPtrB and dump into dXs
  call VecGetOwnershipRange(dXs, istart, iend, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  
  do i=istart, iend-1
     call VecSetValues(dXs, 1, (/i/), (/XsPtrb(i+1)/), INSERT_VALUES, ierr)
     call EChk(ierr, __FILE__, __LINE__)
  end do
  
  ! Must assemble even though we set local values
  call VecAssemblyBegin(dXs, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecAssemblyEnd(dXs, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Deallocate the extra space
  deallocate(XvPtrb, XsPtrb)
  deallocate(mib, Bib, Aib, normalsb, normals0b, xub)

  ! Restore all the arrays
  call VecRestoreArrayF90(XsLocal, XsPtr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecRestoreArrayF90(Xv0, Xv0Ptr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecRestoreArrayF90(Xv, XvPtr, ierr)
  call EChk(ierr,__FILE__,__LINE__)


end subroutine warpDeriv
