subroutine copyRotateVolumeCoordinates()
    use gridData
    implicit none

    ! Working variables
    integer(kind=intType) :: istart, iend, ierr, iRotStart, nDofRot, ii
    integer(kind=intType) :: a1, b1, a2, idxStart
    real(kind=realType) :: Mi(3, 3), coords(3), angle, axis(3)
    real(kind=realType), dimension(:), pointer :: coordsLocal, rotCoords

    ! Get the necessary axisymm info from the grid data
    a1 = mirrorPlaneIdxs(1)  ! Starting index of the mirror plane nodes
    b1 = mirrorPlaneIdxs(2) - 1  ! Ending index of the mirror plane nodes
    a2 = rotPlaneIdxs(1)  ! Starting index of the rotated nodes
    axis = axiSymmAxis  ! Axis of rotation
    angle = axiSymmAngle

    ! Get the local ownership range
    call VecGetOwnershipRange(Xv, istart, iend, ierr)
    call EChk(ierr, __FILE__, __LINE__)
    istart = istart + 1  ! Convert to 1-based indexing
    iend = iend ! Convert to 1-based indexing

    call VecGetArrayF90(Xv, coordsLocal, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    if (istart .ge. a1 .and. istart .lt. b1 .and. iend .le. b1) then
        ! All local volume nodes lie in [a1, b1]
        nDofRot = iend - istart + 1
        idxStart = 1
    else if (istart .ge. a1 .and. istart .lt. b1 .and. iend .gt. b1) then
        ! istart is in [a1, b1] but iend is not
        iend = b1  ! Clip iend to b1
        nDofRot = iend - istart + 1
        idxStart = 1
    else if (istart .lt. a1 .and. iend .gt. a1 .and. iend .lt. b1) then
        ! iend is in [a1, b1] but istart is not
        idxStart = a1 - istart + 1
        istart = a1  ! Clip istart to a1
        nDofRot = iend - istart + 1
    else
        ! The local volume nodes are not in [a1, b1]
        istart = 0
        iend = 0
        nDofRot = 0
    end if

    ! Get the starting index for the rotated coordinates
    iRotStart = a2 + (istart - a1)

    call getRotationMatrixAngleAxis(angle, axis, Mi)

    allocate (rotCoords(nDofRot))

    if (nDofRot .gt. 0) then
        do ii = 1, nDofRot / 3
            ! Extract the x, y, z coordinates
            coords(1) = coordsLocal(idxStart + 3 * ii - 3)
            coords(2) = coordsLocal(idxStart + 3 * ii - 2)
            coords(3) = coordsLocal(idxStart + 3 * ii - 1)

            ! Rodrigues rotation formula
            rotCoords(3 * ii - 2:3 * ii) = matmul(Mi, coords)
        end do

        call VecSetValues(Xv, nDofRot, (/(ii, ii=iRotStart - 1, iRotStart + nDofRot - 2, 1)/), &
                          rotCoords, INSERT_VALUES, ierr)
        call EChk(ierr, __FILE__, __LINE__)
    end if

    ! Need to call VecAssemblyBegin/End
    call VecAssemblyBegin(Xv, ierr)
    call EChk(ierr, __FILE__, __LINE__)
    call VecAssemblyEnd(Xv, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call VecRestoreArrayF90(Xv, coordsLocal, ierr)
    call EChk(ierr, __FILE__, __LINE__)

end subroutine copyRotateVolumeCoordinates

subroutine copyRotateVolumeCoordinates_d()
    use gridData
    implicit none

    ! Working variables
    integer(kind=intType) :: istart, iend, ierr, iRotStart, ndof, ii
    integer(kind=intType) :: a1, b1, a2, idxStart
    real(kind=realType) :: Mi(3, 3), coords_d(3), angle, axis(3)
    real(kind=realType), dimension(:), pointer :: coordsLocal_d, rotCoords_d

    ! Get the necessary axisymm info from the grid data
    a1 = mirrorPlaneIdxs(1)  ! Starting index of the mirror plane nodes
    b1 = mirrorPlaneIdxs(2) - 1  ! Ending index of the mirror plane nodes
    a2 = rotPlaneIdxs(1)  ! Starting index of the rotated nodes
    axis = axiSymmAxis  ! Axis of rotation
    angle = axiSymmAngle

    ! Get the local ownership range
    call VecGetOwnershipRange(dXv, istart, iend, ierr)
    call EChk(ierr, __FILE__, __LINE__)
    istart = istart + 1  ! Convert to 1-based indexing
    iend = iend ! Convert to 1-based indexing

    call VecGetArrayF90(dXv, coordsLocal_d, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    if (istart .ge. a1 .and. istart .lt. b1 .and. iend .le. b1) then
        ! All local volume nodes lie in [a1, b1]
        ndof = iend - istart + 1
        idxStart = 1
    else if (istart .ge. a1 .and. istart .lt. b1 .and. iend .gt. b1) then
        ! istart is in [a1, b1] but iend is not
        iend = b1  ! Clip iend to b1
        nDof = iend - istart + 1
        idxStart = 1
    else if (istart .lt. a1 .and. iend .gt. a1 .and. iend .lt. b1) then
        ! iend is in [a1, b1] but istart is not
        idxStart = a1 - istart + 1
        istart = a1  ! Clip istart to a1
        nDof = iend - istart + 1
    else
        ! The local volume nodes are not in [a1, b1]
        istart = 0
        iend = 0
        ndof = 0
    end if

    ! Get the starting index for the rotated coordinates
    iRotStart = a2 + (istart - a1)

    call getRotationMatrixAngleAxis(angle, axis, Mi)

    allocate (rotCoords_d(ndof))

    if (ndof .gt. 0) then
        do ii = 1, ndof / 3
            ! Extract the x, y, z coordinates
            coords_d(1) = coordsLocal_d(idxStart + 3 * ii - 3)
            coords_d(2) = coordsLocal_d(idxStart + 3 * ii - 2)
            coords_d(3) = coordsLocal_d(idxStart + 3 * ii - 1)

            ! Rodrigues rotation formula
            rotCoords_d(3 * ii - 2:3 * ii) = matmul(Mi, coords_d)
        end do

        call VecSetValues(dXv, ndof, (/(ii, ii=iRotStart - 1, iRotStart + nDof - 2, 1)/), &
                          rotCoords_d, INSERT_VALUES, ierr)
        call EChk(ierr, __FILE__, __LINE__)
    end if

    ! Need to call VecAssemblyBegin/End
    call VecAssemblyBegin(dXv, ierr)
    call EChk(ierr, __FILE__, __LINE__)
    call VecAssemblyEnd(dXv, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call VecRestoreArrayF90(dXv, coordsLocal_d, ierr)
    call EChk(ierr, __FILE__, __LINE__)

end subroutine copyRotateVolumeCoordinates_d

subroutine copyRotateVolumeCoordinates_b()
    use gridData
    implicit none

    ! Working variables
    integer(kind=intType) :: istart, iend, ierr, iRotStart, nDofRot, ii
    integer(kind=intType) :: a1, b1, a2, b2, iStartLocal, totalDof
    real(kind=realType) :: Mi(3, 3), coords_b(3), angle, axis(3)
    real(kind=realType), dimension(:), pointer :: coordsLocal_b
    real(kind=realType), dimension(:), allocatable :: rotCoords_b

    ! Get the necessary axisymm info from the grid data
    a1 = mirrorPlaneIdxs(1)  ! Starting index of the mirror plane nodes
    b1 = mirrorPlaneIdxs(2) - 1  ! Ending index of the mirror plane nodes
    a2 = rotPlaneIdxs(1)  ! Starting index of the rotated nodes
    b2 = rotPlaneIdxs(2) - 1  ! Ending index of the rotated nodes
    axis = axiSymmAxis  ! Axis of rotation
    angle = axiSymmAngle

    ! Get the local ownership range
    call VecGetOwnershipRange(dXv, istart, iend, ierr)
    call EChk(ierr, __FILE__, __LINE__)
    istart = istart + 1  ! Convert to 1-based indexing
    iend = iend ! Convert to 1-based indexing

    call VecGetArrayF90(dXv, coordsLocal_b, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    if (istart .ge. a2 .and. iend .le. b2) then
        ! All local volume nodes lie in [a2, b2]
        nDofRot = iend - istart + 1
        iStartLocal = 1
    else if (istart .lt. a2 .and. iend .gt. a2 .and. iend .le. b2) then
        ! istart is not in [a2, b2] but iend is
        iStartLocal = a2 - istart + 1
        istart = a2  ! Clip istart to a2
        nDofRot = iend - istart + 1
    else if (istart .ge. a2 .and. istart .lt. b2 .and. iend .gt. b2) then
        ! istart is in [a2, b2] but iend is not
        iend = b2  ! Clip iend to b2
        nDofRot = iend - istart + 1
        iStartLocal = 1
    else
        istart = 0
        iend = 0
        nDofRot = 0
    end if

    ! Get the starting index for the rotated coordinates
    iRotStart = a1 + (istart - a2)
    allocate (rotCoords_b(nDofRot))

    call getRotationMatrixAngleAxis(-angle, axis, Mi)

    if (nDofRot .gt. 0) then
        do ii = 1, nDofRot / 3
            ! Extract the x, y, z coordinates
            coords_b(1) = coordsLocal_b(iStartLocal + 3 * ii - 3)
            coords_b(2) = coordsLocal_b(iStartLocal + 3 * ii - 2)
            coords_b(3) = coordsLocal_b(iStartLocal + 3 * ii - 1)

            ! Rodrigues rotation formula
            ! Add the rotated seeds to the original seeds
            rotCoords_b(3 * ii - 2:3 * ii) = matmul(Mi, coords_b)
        end do

        ! Set the rotated + original seeds in the correct place
        call VecSetValues(dXv, nDofRot, (/(ii, ii=iRotStart - 1, iRotStart + nDofRot - 2)/), &
                          rotCoords_b, ADD_VALUES, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! Zero out all seeds on the rotated plane
        coordsLocal_b(iStartLocal:iStartLocal + nDofRot - 1) = zero
    end if

    ! Need to call VecAssemblyBegin/End
    call VecAssemblyBegin(dXv, ierr)
    call EChk(ierr, __FILE__, __LINE__)
    call VecAssemblyEnd(dXv, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Restore the array to deallocate the coordsLocal_b pointer
    call VecRestoreArrayF90(dXv, coordsLocal_b, ierr)
    call EChk(ierr, __FILE__, __LINE__)

end subroutine copyRotateVolumeCoordinates_b
