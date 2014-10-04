subroutine warpMesh()

  ! This is the master mesh warping gateway function. It does not do
  ! any actual computations, rather just extracts the required
  ! pointers from petsc, determines which of the two gateway functions
  ! are called depending on if the fast or exact code are called.

  use gridData
  use gridInput
  use communication
  implicit none

  ! Working parameters
  integer(kind=intType) :: ierr

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

  ! Now call required routine:
  if (evalMode == EVAL_EXACT) then 
     call warpMeshExact()
  else if (evalMode == EVAL_FAST) then
     call warpMeshFast()
  end if
     
  ! Restore all the arrays
  call VecRestoreArrayF90(XsLocal, XsPtr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecRestoreArrayF90(Xv0, Xv0Ptr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecRestoreArrayF90(Xv, XvPtr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

end subroutine warpMesh

