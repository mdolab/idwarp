subroutine warpDeriv(dXv_f, ndof_warp)

  use gridData
  use gridInput
  use communication
  implicit none

  ! Input
  integer(kind=intType), intent(in) :: ndof_warp
  real(kind=realType), intent(in) :: dXv_f(ndof_warp)

  ! ! Working Data
  ! integer(kind=intType) :: i, istart, iend, ierr, isize
  ! real(kind=realType), dimension(:), allocatable :: dxsSummed
  ! ! Scatter Xs into our local vector  
  ! call VecScatterBegin(XsToXsLocal, Xs, XsLocal, INSERT_VALUES, SCATTER_FORWARD, ierr)
  ! call EChk(ierr,__FILE__,__LINE__)

  ! call VecScatterEnd(XsToXsLocal, Xs, XsLocal, INSERT_VALUES, SCATTER_FORWARD, ierr)
  ! call EChk(ierr,__FILE__,__LINE__)

  ! ! Extract a pointer from Xs to use in the main routine
  ! call VecGetArrayF90(XsLocal, XsPtr, ierr)
  ! call EChk(ierr,__FILE__,__LINE__)

  ! ! Extract the pointers for the volume nodes that we wish to operate on
  ! call VecGetArrayF90(Xv0, Xv0Ptr, ierr)
  ! call EChk(ierr,__FILE__,__LINE__)

  ! call VecGetArrayF90(Xv, XvPtr, ierr)
  ! call EChk(ierr,__FILE__,__LINE__)

  ! call VecGetArrayF90(dXs, dXsPtr, ierr)
  ! call EChk(ierr,__FILE__,__LINE__)

  ! ! Allocate the extra data we need for the warping derivative:
  ! allocate(XvPtrb(size(XvPtr)), XsPtrb(SIZE(XsPtr)))
  ! allocate(xub(3, nUnique))
  ! allocate(Mib(3, 3, nUnique), Bib(3, nUnique))
  ! allocate(normalsb(3, nUnique), normals0b(3, nUnique))
 
  ! ! dXv_f is the actual reverse seed so copy:
  ! XvPtrb(:) = dXv_f

  ! ! Zero the output surface derivative vector dXs
  ! call vecZeroEntries(dXs, ierr)
  ! call EChk(ierr, __FILE__, __LINE__)

  ! ! Now run the actual reverse mode routine
  ! if (evalMode == EVAL_EXACT) then 
  !    call warpMeshExact_b()
  ! else
  !    call warpMeshFast_b()
  ! end if

  ! call VecGetSize(dXs, isize, ierr)
  ! call EChk(ierr, __FILE__, __LINE__)

  ! call VecGetOwnershipRange(Xs, istart, iend, ierr)
  ! call EChk(ierr, __FILE__, __LINE__)

  ! allocate(dxssummed(isize))
  ! dxssummed = zero

  ! call MPI_Allreduce(xsptrb, dxssummed, isize, MPI_REAL8, MPI_SUM, warp_comm_world, ierr)
  ! call EChk(ierr, __FILE__, __LINE__)
  
  ! ! Copy just the required part into dXsPtr
  ! dXsPtr = dxssummed(iStart+1:iEnd)
  ! deallocate(dxssummed)
  
  ! ! Deallocate the extra space
  ! deallocate(XvPtrb, XsPtrb)
  ! deallocate(mib, Bib, normalsb, normals0b, xub)

  ! ! Restore all the arrays
  ! call VecRestoreArrayF90(XsLocal, XsPtr, ierr)
  ! call EChk(ierr,__FILE__,__LINE__)

  ! call VecRestoreArrayF90(Xv0, Xv0Ptr, ierr)
  ! call EChk(ierr,__FILE__,__LINE__)

  ! call VecRestoreArrayF90(Xv, XvPtr, ierr)
  ! call EChk(ierr,__FILE__,__LINE__)

  ! call VecRestoreArrayF90(dXs, dXsPtr, ierr)
  ! call EChk(ierr,__FILE__,__LINE__)
  
end subroutine warpDeriv
