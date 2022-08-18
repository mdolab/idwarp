subroutine initPETSc(comm)
    ! Simple wrapper for initializing petsc in case it wasn't done in
    ! python. This also gets the comm data.
    use communication
    use gridData
    implicit none

    ! Input params
    integer(kind=intType) :: comm

    ! Working variables
    integer(kind=intType) :: ierr

    call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    warp_comm_world = comm
    warp_comm_self = mpi_comm_self
    call MPI_Comm_size(warp_comm_world, nProc, ierr)
    call MPI_Comm_rank(warp_comm_world, myid, ierr)

end subroutine initPETSc

