subroutine copyRotateVolumeCoordinates()
    use gridData
    implicit none

    ! Working variables
    integer(kind=intType) :: istart, iend, ierr, iRotStart, ndof, ii
    integer(kind=intType) :: a1, b1, a2
    real(kind=realType) :: Mi(3, 3), coords(3), angle, axis(3)
    real(kind=realType), dimension(:), pointer :: temp, rotCoords
    integer(kind=intType), dimension(:), allocatable :: idx

    ! Get the necessary axisymm info from the grid data
    a1 = mirrorPlaneIdxs(1)  ! Starting index of the mirror plane nodes
    b1 = mirrorPlaneIdxs(2)  ! Ending index of the mirror plane nodes
    a2 = rotPlaneIdxs(1)  ! Starting index of the rotated nodes
    axis = axiSymmAxis  ! Axis of rotation
    angle = axiSymmAngle

    ! Get the local ownership range
    call VecGetOwnershipRange(Xv, istart, iend, ierr)
    call EChk(ierr, __FILE__, __LINE__)
    istart = istart + 1  ! Convert to 1-based indexing
    iend = iend + 1 ! Convert to 1-based indexing

    call VecGetArrayF90(Xv, temp, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    if (istart .ge. a1 .and. istart .le. b1) then
        ! All local volume nodes lie in [a1, b1]
        ndof = size(temp)
    else if (istart .gt. a1 .and. istart .lt. iend .and. iend .gt. b1) then
        ! istart is in [a1, b1] but iend is not
        iend = b1  ! Clip the end to b1
        ndof = size(temp(istart:iend))
    else if (istart .le. a1 .and. iend .gt. a1 .and. iend .lt. b1) then
        ! iend is in [a1, b1] but istart is not
        istart = a1  ! Clip the start to a1
        ndof = size(temp(istart:iend))
    else
        ! The local volume nodes are not in [a1, b1]
        istart = 0
        iend = 0
        ndof = 0
    end if

    ! Get the starting index for the rotated coordinates
    iRotStart = a2 + (istart - a1)

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
        idx = idx + iRotStart

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
    use gridData
    implicit none

    ! Working variables
    integer(kind=intType) :: istart, iend, ierr, iRotStart, ndof, ii
    integer(kind=intType) :: a1, b1, a2
    real(kind=realType) :: Mi(3, 3), coords(3), angle, axis(3)
    real(kind=realType), dimension(:), pointer :: temp, rotCoords
    integer(kind=intType), dimension(:), allocatable :: idx

    ! Get the necessary axisymm info from the grid data
    a1 = mirrorPlaneIdxs(1)  ! Starting index of the mirror plane nodes
    b1 = mirrorPlaneIdxs(2)  ! Ending index of the mirror plane nodes
    a2 = rotPlaneIdxs(1)  ! Starting index of the rotated nodes
    axis = axiSymmAxis  ! Axis of rotation
    angle = axiSymmAngle

    ! Get the local ownership range
    call VecGetOwnershipRange(dXv, istart, iend, ierr)
    call EChk(ierr, __FILE__, __LINE__)
    istart = istart + 1  ! Convert to 1-based indexing
    iend = iend + 1 ! Convert to 1-based indexing

    call VecGetArrayF90(dXv, temp, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    if (istart .ge. a1 .and. istart .le. b1) then
        ! All local volume nodes lie in [a1, b1]
        ndof = size(temp)
    else if (istart .gt. a1 .and. istart .lt. iend .and. iend .gt. b1) then
        ! istart is in [a1, b1] but iend is not
        iend = b1  ! Clip the end to b1
        ndof = size(temp(istart:iend))
    else if (istart .le. a1 .and. iend .gt. a1 .and. iend .lt. b1) then
        ! iend is in [a1, b1] but istart is not
        istart = a1  ! Clip the start to a1
        ndof = size(temp(istart:iend))
    else
        ! The local volume nodes are not in [a1, b1]
        istart = 0
        iend = 0
        ndof = 0
    end if

    ! Get the starting index for the rotated coordinates
    iRotStart = a2 + (istart - a1)

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

    ! Move the indices to the correct location in the full dXv array
    if (ndof .gt. 0) then
        idx = idx - 2
        idx = idx + iRotStart

        call VecSetValues(dXv, ndof, idx, rotCoords, INSERT_VALUES, ierr)
        call EChk(ierr, __FILE__, __LINE__)
    end if

    ! Need to call VecAssemblyBegin/End
    call VecAssemblyBegin(dXv, ierr)
    call EChk(ierr, __FILE__, __LINE__)
    call VecAssemblyEnd(dXv, ierr)
    call EChk(ierr, __FILE__, __LINE__)

end subroutine copyRotateVolumeSeeds

! Need to get a mapping of the load balancing (rank, istart, iend)
! Pass the owners