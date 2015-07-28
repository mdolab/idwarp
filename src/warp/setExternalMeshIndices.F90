subroutine setExternalMeshIndices(ndof_solver, solver_indices)

  use communication
  use gridData
  implicit none

  ! Subroutine Arguments
  integer(kind=intType), intent(in)  :: ndof_solver
  integer(kind=intType), intent(in)  :: solver_indices(ndof_solver)

  ! Local Variables
  integer(kind=intType) :: cumDOFSolverProc(0:nProc)
  integer(kind=intType) :: DOFSolverProc(nProc)
  integer(kind=intType) :: ierr, iproc

  ! Set the number of DOF for the solver
  solvermeshdof = ndof_solver

  ! Compute the cumDOFProc for the solver
  call mpi_allgather(ndof_solver, 1, MPI_INTEGER, &
       DOFSolverProc, 1, mpi_integer4, warp_comm_world, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  
  cumDOFSolverProc(0) = 0
  do iproc=1, nproc
     cumDOFSolverProc(iproc) = cumDOFSolverProc(iproc-1) + DOFSolverProc(iproc)
  end do

  ! Create an 'empty' PETSc vector. This will have an arrary
  ! vecPlace'd in them for doing the actual scatters.
  call VecCreateMPIWithArray(WARP_COMM_WORLD, 1, ndof_solver, PETSC_DECIDE, &
       PETSC_NULL_SCALAR, solverGridVec, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! ~~~~~~~~ Common to Solver Scatter~~~~~~~~~
  call ISCreateGeneral(warp_comm_world, ndof_solver, solver_indices, &
       PETSC_COPY_VALUES, IS1, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  call ISCreateStride (warp_comm_world, ndof_solver, &
       cumDOFSolverProc(myid), 1, IS2, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecScatterCreate(commonGridVec, IS1, solverGridVec, IS2, common_to_solver, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  
  call ISDestroy(IS1, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call ISDestroy(IS2, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  gridIndicesSet = 1

end subroutine setExternalMeshIndices
