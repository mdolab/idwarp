subroutine initNodalProperties

  ! This wrapper routine is to be called from python after the surface
  ! mesh is set to compute the initial nodal properties. 

  use gridData
  use gridInput
  use communication
  implicit none

  integer(kind=intType) :: ierr

  ! Scatter Xs into our local vector  
  call VecScatterBegin(XsToXsLocal, Xs, XsLocal, INSERT_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecScatterEnd(XsToXsLocal, Xs, XsLocal, INSERT_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Extract a pointer from Xs to use in the main routine
  call VecGetArrayF90(XsLocal, XsPtr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call computeNodalProperties(.True.)

  ! Restore all the arrays
  call VecRestoreArrayF90(XsLocal, XsPtr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

end subroutine initNodalProperties
