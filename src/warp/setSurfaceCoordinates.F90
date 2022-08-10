subroutine setSurfaceCoordinates(coordinates, cdof)

    use gridData
    implicit none

    ! Input Arguments
    integer(kind=intType), intent(in) :: cdof
    real(kind=realType), intent(in) :: coordinates(cdof)

    ! Local Arguments
    integer(kind=intType) :: ierr, istart, iend, i

    call VecGetOwnershipRange(Xs, istart, iend, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    do i = 1, cdof
        call VecSetValues(Xs, 1, [iStart + i - 1], [coordinates(i)], INSERT_VALUES, ierr)
        call EChk(ierr, __FILE__, __LINE__)
    end do

    ! While we only set local values, we STILL have to call
    ! VecAssemblyBegin/End
    call VecAssemblyBegin(Xs, ierr)
    call EChk(ierr, __FILE__, __LINE__)
    call VecAssemblyEnd(Xs, ierr)
    call EChk(ierr, __FILE__, __LINE__)

end subroutine setSurfaceCoordinates
