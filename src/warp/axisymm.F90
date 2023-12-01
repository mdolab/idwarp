subroutine copyRotateVolumeCoordinates(iMirrorStart, angle, axis)
    use gridData
    use communication
    implicit none

    ! Input variables
    integer(kind=intType), intent(in) :: iMirrorStart
    real(kind=realType), intent(in) :: angle, axis(3)

    ! Working variables
    integer(kind=intType) :: istart, iend, ierr, iRotStart, ndof, ii, nSizeGlobal, iHalf
    real(kind=realType) :: Mi(3, 3), coords(3)
    real(kind=realType), dimension(:), pointer :: temp, rotCoords
    integer(kind=intType), dimension(:), allocatable :: idx

    ! Get the local ownership range
    call VecGetOwnershipRange(Xv, istart, iend, ierr)
    call EChk(ierr, __FILE__, __LINE__)
    istart = istart + 1  ! Need to increment to one-based indexing

    call VecGetSize(Xv, nSizeGlobal, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    iHalf = nSizeGlobal / 2 ! integer division

    call VecGetArrayF90(Xv, temp, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Check is the entire vector is on the mirror plane
    if (istart .lt. iHalf .and. iend .le. iHalf) then
        ndof = size(temp)
    else if (istart .le. iHalf .and. iend .gt. iHalf) then
        ! Vector is partially on mirror and rotated planes
        iend = iHalf  ! truncate to only work with the mirror plane
        ndof = iend - istart + 1
    else
        ndof = 0
    end if

    call getRotationMatrixAngleAxis(angle, axis, Mi)

    allocate (idx(ndof))
    allocate (rotCoords(ndof))

    do ii = 1, ndof / 3
        ! Extract the x, y, z coordinates
        coords(1) = temp(3 * ii - 2)
        coords(2) = temp(3 * ii - 1)
        coords(3) = temp(3 * ii)

        ! Rodrigues rotation formula
        rotCoords(3 * ii - 2:3 * ii) = matmul(Mi, coords)

        ! Set the indices of the rotated coordinates
        idx(3 * ii - 2) = 3 * ii - 2
        idx(3 * ii - 1) = 3 * ii - 1
        idx(3 * ii) = 3 * ii
    end do

    ! Move the indices to the correct location in the full Xv array
    if (ndof .gt. 0) then
        idx = idx - 2
        idx = idx + iStart + iHalf

        call VecSetValues(Xv, ndof, idx, rotCoords, INSERT_VALUES, ierr)
        call EChk(ierr, __FILE__, __LINE__)
    end if

    ! Need to call VecAssemblyBegin/End
    call VecAssemblyBegin(Xv, ierr)
    call EChk(ierr, __FILE__, __LINE__)
    call VecAssemblyEnd(Xv, ierr)
    call EChk(ierr, __FILE__, __LINE__)

end subroutine copyRotateVolumeCoordinates

subroutine copyRotateVolumeSeeds()
end subroutine copyRotateVolumeSeeds
