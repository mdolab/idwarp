subroutine createCommonGrid(volNodes, nVolLocal)

    ! This routine create a few arrys that are the size of the "common"
    ! grid...that is the ordering that is given in the original grid
    ! file. Generally this will result in a dumb "nVol/nProc"
    ! distribution, while exact, does not result in proper load
    ! balancing since each volume node will take a different amount of
    ! time.

    use gridData
    use communication
    implicit none

    ! Subroutine variables
    integer(kind=intType), intent(in) :: nVolLocal
    real(kind=realType), dimension(3, nVolLocal), intent(in) :: volNodes

    ! Working variables
    integer(kind=intType) :: i, j, nVol, ierr, nSurface
    integer(kind=intType), dimension(:), allocatable :: volNodesProc

    ! Gather the displacements
    allocate (volNodesProc(0:nProc))
    volNodesProc(:) = 0
    call MPI_allgather(nVolLocal, 1, MPI_INTEGER, volNodesProc(1:), 1, MPI_INTEGER, &
                       warp_comm_world, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Finish the displacment calc:
    do i = 2, nProc
        volNodesProc(i) = volNodesProc(i) + volNodesProc(i - 1)
    end do
    nVol = volNodesProc(nProc)

    if (myid == 0) then
        write (*, "(a)", advance="no") '#------------------------------#'
        print "(1x)"
        write (*, "(a)", advance="no") " Total Volume Nodes : "
        write (*, "(I9,1x)", advance="no") nVol
        print "(1x)"
        write (*, "(a)", advance="no") '#------------------------------#'
        print "(1x)"
    end if

    ! Now create the master Xv array:
    call VecCreate(warp_comm_world, commonGridVec, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Set to be be blocked
    call VecSetBlockSize(commonGridVec, 3, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Type and size
    call VecSetType(commonGridVec, "mpi", ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call VecSetSizes(commonGridVec, nVolLocal * 3, PETSC_DECIDE, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Now each processor adds it's own nodes (It is possible only one proc does it)
    do i = 1, nVolLocal
        call VecSetValuesBlocked(commonGridVec, 1, (/volNodesProc(myid) + i - 1/), volNodes(:, i), INSERT_VALUES, ierr)
        call EChk(ierr, __FILE__, __LINE__)
    end do

    call VecAssemblyBegin(commonGridVec, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call VecAssemblyEnd(commonGridVec, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    deallocate (volNodesProc)
    commonGridVecSet = 1

    commonMeshDOF = nVol * 3

end subroutine createCommonGrid
