! This file implements the fake mesh warp -- it only takes peturbation
! of the degrees of freedom that are already on the surface. 
subroutine warpDerivSurfOnly(dXv_f, ndof)

  use gridData
  implicit none

  ! Input
  integer(kind=intType), intent(in) :: ndof
  real(kind=realType), intent(in) :: dXv_f(ndof)

  ! Local Variables
  integer(kind=intType) :: ierr

  ! Place arrays
  call VecPlaceArray(solverGridVec, dXv_f, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  
  call VecZeroEntries(commonGridVec, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Now do the vecScatter: common_to_solver in REVERSE
  call VecScatterBegin(common_to_solver, solverGridVec, commonGridVec, &
       ADD_VALUES, SCATTER_REVERSE, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  call VecScatterEnd  (common_to_solver, solverGridVec, commonGridVec, &
       ADD_VALUES, SCATTER_REVERSE, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! And now scatter from commonGridVec to dXs in FORWARD
  call VecZeroEntries(dXs, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecScatterBegin(common_to_dXs, commonGridVec, dXs, &
       ADD_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  call VecScatterEnd  (common_to_dXs, solverGridVec, dXs, &
       ADD_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  
  call VecResetArray(solverGridVec, ierr)
  call EChk(ierr, __FILE__, __LINE__)
end subroutine warpDerivSurfOnly
