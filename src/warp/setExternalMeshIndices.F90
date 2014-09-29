subroutine setExternalMeshIndices(ndof_solver, solver_indices)
  use communication
  use gridData
  implicit none
#include "include/petscversion.h"
  ! Subroutine Arguments
  integer(kind=intType), intent(in)  :: ndof_solver
  integer(kind=intType), intent(in)  :: solver_indices(ndof_solver)

  ! Local Variables
  integer(kind=intType) :: iLow, iHigh
  integer(kind=intType) :: cumDOFSolverProc(0:nProc)
  integer(kind=intType) :: DOFSolverProc(nProc)
  integer(kind=intType) :: ierr, iproc
  integer(kind=intType) :: ndof_warp, totalDOF

  ! Compute the cumDOFProc for the solver
  call mpi_allgather(ndof_solver, 1, MPI_INTEGER, &
       DOFSolverProc, 1, mpi_integer4, warp_comm_world, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  cumDOFSolverProc(0) = 0
  do iproc=1, nproc
     cumDOFSolverProc(iproc) = cumDOFSolverProc(iproc-1) + DOFSolverProc(iproc)
  end do

  ! Get the number of owned nodes here:
  call VecGetOwnershipRange(Xv, iLow, iHigh, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  ndof_warp = iHigh-iLow

  ! Total mesh size
  call VecGetSize(Xv, totalDOF, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  
  ! Setup the PETSc Vectors
  ! cgnsGrid Vector
  if (.not. commonGridVecSet) then
     call VecCreate(warp_comm_world, commonGridVec, ierr)
     call EChk(ierr, __FILE__, __LINE__)
     call VecSetType(commonGridVec, "mpi", ierr)
     call EChk(ierr, __FILE__, __LINE__)
     call VecSetSizes(commonGridVec, PETSC_DECIDE, totalDOF, ierr)
     call EChk(ierr, __FILE__, __LINE__)
  end if

  ! Create two 'empty' PETSc vectors. These will have an arrary
  ! vecPlace'd in them for doing the actual scatters. 
  call VecCreateMPIWithArray(WARP_COMM_WORLD, 1, ndof_solver, PETSC_DECIDE, &
       PETSC_NULL_SCALAR, solverGridVec, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  
  call VecCreateMPIWithArray(WARP_COMM_WORLD, 1, ndof_warp, PETSC_DECIDE, &
       PETSC_NULL_SCALAR, dXv, ierr)
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

  ! ~~~~~~~~ Common to Warp Scatter ~~~~~~~~~~~~

  ! Note that currently, the warp grid is identical to the CGNS
  ! grid. This may be modified at a later time. 
  call ISCreateStride(warp_comm_world, ndof_warp, iLow, 1, IS1, ierr)
  call EChk(ierr, __FILE__, __LINE__)  

  call ISCreateStride (warp_comm_world, ndof_warp, iLow, 1, IS2, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecScatterCreate(commonGridVec, IS1, dXv, IS2, common_to_warp, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call ISDestroy(IS1, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  
  call ISDestroy(IS2, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  gridIndicesSet = .True.
  commonGridVecSet = .True.

end subroutine setExternalMeshIndices
