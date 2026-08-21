! This routine get/sets the current grid nodes in Xv to python
subroutine getVolumeCoordinates(gridNodes, nDOF)

    use gridData
    use petscCompat, only: VecGetArrayCompat, VecRestoreArrayCompat
    implicit none

    ! Subroutine variables
    integer(kind=intType), intent(in) :: nDOF
    real(kind=realType), dimension(nDOF), intent(out) :: gridNodes
    real(kind=realType), pointer, dimension(:) :: xx
    integer(kind=intType) :: ierr

    ! No error checking here!
    call VecGetArrayCompat(Xv, xx, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Perform copy
    gridNodes = xx

    call VecRestoreArrayCompat(Xv, xx, ierr)
    call EChk(ierr, __FILE__, __LINE__)

end subroutine getVolumeCoordinates

subroutine setVolumeCoordinates(gridNodes, nDOF)

    use gridData
    use petscCompat, only: VecGetArrayCompat, VecRestoreArrayCompat
    implicit none

    ! Subroutine variables
    integer(kind=intType), intent(in) :: nDOF
    real(kind=realType), dimension(nDOF), intent(in) :: gridNodes
    real(kind=realType), pointer, dimension(:) :: xx
    integer(kind=intType) :: ierr

    ! No error checking here!
    call VecGetArrayCompat(Xv, xx, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Perform copy
    xx = gridNodes

    call VecRestoreArrayCompat(Xv, xx, ierr)
    call EChk(ierr, __FILE__, __LINE__)

end subroutine setVolumeCoordinates
