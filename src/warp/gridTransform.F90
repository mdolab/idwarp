! This file contains two routines that transform the "warp" grid to
! the "solver" grid and vice versa

subroutine warp_to_solver_grid(warp_grid, wdof, solver_grid, sdof)

  use gridData
  implicit none

  ! Subroutine arguments
  integer(kind=intType), intent(in)  :: wdof, sdof
  real(kind=realType), intent(in)   :: warp_grid(wdof)
  real(kind=realType), intent(inout)   :: solver_grid(sdof)

  ! Local Variables
  integer(kind=intType) :: ierr

  ! Place Arrays in vectors
  call VecPlaceArray(dXv, warp_grid, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecPlaceArray(solverGridVec, solver_grid, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  
  call VecZeroEntries(commonGridVec, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Now do the vecScatter: common_to_warp in REVERSE
  call VecScatterBegin(common_to_warp, dXv, commonGridVec, &
       INSERT_VALUES, SCATTER_REVERSE, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecScatterEnd  (common_to_warp, dXv, commonGridVec, &
       INSERT_VALUES, SCATTER_REVERSE, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Now do the vecScatter: common_to_solver in FORWARD
  call VecScatterBegin(common_to_solver, commonGridVec, solverGridVec, &
       INSERT_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  call VecScatterEnd  (common_to_solver, commonGridVec, solverGridVec, &
       INSERT_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Destroy the vectors
  call VecResetArray(dXv, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecResetArray(solverGridVec, ierr)
  call EChk(ierr, __FILE__, __LINE__)

end subroutine warp_to_solver_grid

subroutine solver_to_warp_grid(solver_grid, sdof, warp_grid, wdof)

  use gridData
  implicit none

  ! Subroutine arguments
  integer(kind=intType), intent(in)  :: wdof, sdof
  real(kind=realType), intent(in)   :: solver_grid(sdof)
  real(kind=realType), intent(inout)   :: warp_grid(wdof)

  ! Local Variables
  integer(kind=intType) :: ierr

  ! Place arrays
  call VecPlaceArray(solverGridVec, solver_grid, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  
  call VecPlaceArray(dXv, warp_grid, ierr)
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

  ! Now do the vecScatter: common_to_warp in FOREWARD
  call VecScatterBegin(common_to_warp, commonGridVec, dXv, &
       INSERT_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  call VecScatterEnd  (common_to_warp, commonGridVec, dXv, &
       INSERT_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Reset Arrays
  call VecResetArray(dXv, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecResetArray(solverGridVec, ierr)
  call EChk(ierr, __FILE__, __LINE__)

end subroutine solver_to_warp_grid
