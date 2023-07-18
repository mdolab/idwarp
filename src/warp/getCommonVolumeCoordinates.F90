! This routine get/sets the current grid nodes from the Common grid
! ordering.
subroutine getCommonVolumeCoordinates(gridNodes, nDOF)

    use gridData
    implicit none

    ! Subroutine variables
    integer(kind=intType), intent(in) :: nDOF
    real(kind=realType), dimension(nDOF), intent(out) :: gridNodes
    real(kind=realType), pointer, dimension(:) :: xx
    integer(kind=intType) :: ierr

    ! Now do the vecScatter: warp_to_common in REVERSE
    call VecScatterBegin(common_to_warp, Xv, commonGridVec, &
                         INSERT_VALUES, SCATTER_REVERSE, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call VecScatterEnd(common_to_warp, Xv, commonGridVec, &
                       INSERT_VALUES, SCATTER_REVERSE, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call VecGetArrayF90(commonGridVec, xx, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Perform the actual copy
    gridNodes = xx

    call VecRestoreArrayF90(Xv, xx, ierr)
    call EChk(ierr, __FILE__, __LINE__)

end subroutine getCommonVolumeCoordinates

subroutine setCommonVolumeCoordinates(gridNodes, nDOF)

    use gridData
    implicit none

    ! Subroutine variables
    integer(kind=intType), intent(in) :: nDOF
    real(kind=realType), dimension(nDOF), intent(in) :: gridNodes
    real(kind=realType), pointer, dimension(:) :: xx
    integer(kind=intType) :: ierr

    call VecGetArrayF90(commonGridVec, xx, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Perform actual copy
    xx = gridNodes

    call VecRestoreArrayF90(commonGridVec, xx, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Now do the vecScatter: warp_to_common in REVERSE
    call VecScatterBegin(common_to_warp, commonGridVec, Xv, &
                         INSERT_VALUES, SCATTER_FORWARD, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call VecScatterEnd(common_to_warp, commonGridVec, Xv, &
                       INSERT_VALUES, SCATTER_FORWARD, ierr)
    call EChk(ierr, __FILE__, __LINE__)

end subroutine setCommonVolumeCoordinates

