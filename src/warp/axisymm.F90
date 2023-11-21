subroutine copyRotateVolumeCoordinates(iMirrorStart, iMirrorEnd, angle, axis)
    use gridData
    use communication
    implicit none

    ! Input variables
    integer(kind=intType), intent(in) :: iMirrorStart, iMirrorEnd
    real(kind=realType), intent(in) :: angle, axis(3)

    ! Working variables
    integer(kind=intType) :: istart, iend, ierr, iRotStart, ndof, ii, nVol
    real(kind=realType) :: x, y, z, a, b, c, overNorm, sin_theta, cos_theta
    real(kind=realType), dimension(:), pointer :: temp
    integer(kind=intType), dimension(:), allocatable :: idx

    ! Get the local ownership range
    call VecGetOwnershipRange(Xv, istart, iend, ierr)
    call EChk(ierr, __FILE__, __LINE__)
    istart = istart + 1
    iend = iend + 1

    call VecGetArrayF90(Xv, temp, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    nVol = size(temp) / 3

    a = axis(1)
    b = axis(2)
    c = axis(3)

    overNorm = 1.0 / sqrt(a**2 + b**2 + c**2)
    a = a * overNorm
    b = b * overNorm
    c = c * overNorm

    cos_theta = cos(angle)
    sin_theta = sin(angle)

    ! Need to determine if we have mirror plane points
    ! Case 1: All local points lie on the mirror plane
    ndof = 0
    iRotStart = 0
    if (iMirrorStart < iend .and. iend .le. iMirrorEnd) then
        print *, myID, "Loop 1"
        ndof = iend - istart
        iRotStart = iMirrorEnd + istart - 1
    end if

    ! Case 2: Some mirror points lie on the mirror plane
    if (iMirrorStart < iend .and. iMirrorEnd < iend) then
        print *, myID, "Loop 2"
        ndof = iMirrorEnd - istart
        iRotStart = iMirrorEnd + istart - 1
    end if

    if (ndof .eq. 0 .and. iRotStart .eq. 0) then
        return
    end if

    print *, myID, istart, iend, iend - istart
    print *, myID, iMirrorStart, iMirrorEnd
    print *, myID, ndof, iRotStart

    allocate (idx(ndof))

    do ii = 1, nVol
        ! Extract the x, y, z coordinates
        x = temp(3 * ii - 2)
        y = temp(3 * ii - 1)
        z = temp(3 * ii)

        temp(3 * ii - 2) = a**2 * (1 - cos_theta) * x + &
                           (b * a * (1 - cos_theta) - c * sin_theta) * y + &
                           (c * a * (1 - cos_theta) + b * sin_theta) * z

        temp(3 * ii - 1) = (a * b * (1 - cos_theta) + c * sin_theta) * x + &
                           b**2 * (1 - cos_theta) * y + &
                           (c * b * (1 - cos_theta) - a * sin_theta) * z

        temp(3 * ii) = (a * c * (1 - cos_theta) - b * sin_theta) * x + &
                       (b * c * (1 - cos_theta) + a * sin_theta) * y + &
                       c**2 * (1 - cos_theta) * z

        idx(3 * ii - 2) = 3 * ii - 2
        idx(3 * ii - 1) = 3 * ii - 1
        idx(3 * ii) = 3 * ii
    end do

    idx = iRotStart + idx - 1
    print *, idx(1), idx(size(idx)), idx(size(idx)) - idx(1)

    call VecSetValues(Xv, ndof, idx, temp(1:ndof), INSERT_VALUES, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Need to call VecAssemblyBegin/End
    call VecAssemblyBegin(Xv, ierr)
    call EChk(ierr, __FILE__, __LINE__)
    call VecAssemblyEnd(Xv, ierr)
    call EChk(ierr, __FILE__, __LINE__)
end subroutine copyRotateVolumeCoordinates
