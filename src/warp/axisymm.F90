subroutine copyRotateVolumeCoordinates()
    use communication
    use precision
    use gridData
    use cgnsGrid
    implicit none

    ! Working variables
    integer(kind=intType) :: i, j, k, ii, jj, istart, iend, ierr, iiLocal
    integer(kind=intType) :: iZone, nZones, iBoco, ptRange(3, 2), rotDir, iPlane
    integer(kind=intType) :: iNode, dir1, dir2, rotOffset, nNodes, mirrorIdx
    integer(kind=intType), dimension(:), allocatable :: insertIdxs
    real(kind=realType) :: Mi(3, 3), coords(3), rotCoords(3)
    real(kind=realType), dimension(:), allocatable :: insertCoords
    real(kind=realType), dimension(:), pointer :: coordsLocal

    nZones = size(zones)

    ! Get the rotation matrix
    call getRotationMatrixAngleAxis(axiSymmAngle, axiSymmAxis, Mi)

    ! Get the range of the volume nodes on this proc
    call VecGetOwnershipRange(Xv, istart, iend, ierr)
    istart = istart + 1  ! Convert to 1-based indexing
    call EChk(ierr, __FILE__, __LINE__)

    ! Get the local array of coordinates
    call VecGetArrayF90(Xv, coordsLocal, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ii = 0
    jj = 0
    ! Loop over the zones in the grid
    zoneLoop: do iZone = 1, nZones
        rotDir = 0 ! Reset the rotation direction

        ! Loop over the BCs to see if there is a mirror plane BC
        bocoLoop: do iBoco = 1, size(zones(iZone)%bocos)
            ! Check if the BCs for this zone contain the mirror family
            if (zones(iZone)%bocos(iBoco)%family .eq. mirrorFamily) then
                ! Determine the rotation direction
                ! This is the direction in which the zone is one cell wide
                ! The point range for the BC will be equal in this direction
                ! signifying that the zone is one cell wide in this direction
                ptRange = zones(iZone)%bocos(iBoco)%ptRange
                if (ptRange(1, 1) == ptRange(1, 2)) then
                    mirrorIdx = ptRange(1, 1)
                    rotDir = 1 ! i-direction
                else if (ptRange(2, 1) == ptRange(2, 2)) then
                    mirrorIdx = ptRange(2, 1)
                    rotDir = 2 ! j-direction
                else if (ptRange(3, 1) == ptRange(3, 2)) then
                    mirrorIdx = ptRange(3, 1)
                    rotDir = 3 ! k-direction
                end if
            end if
        end do bocoLoop

        if (rotDir == 0) then
            ! No mirror plane BC found, skip this zone
            cycle zoneLoop
        end if

        ! Determine the plane directions to loop over and set the rotation
        ! plane offset (These come from the loop structure in readStructuredCGNS.F90)
        if (rotDir == 1) then
            ! The zone is one cell wide in the i-direction (final loop in readStructuredCGNS.F90)
            dir1 = zones(izone)%jl
            dir2 = zones(izone)%kl
            rotOffset = 1
        else if (rotDir == 2) then
            ! The zone is one cell wide in the j-direction (second loop in readStructuredCGNS.F90)
            dir1 = zones(iZone)%il
            dir2 = zones(iZone)%kl
            rotOffset = zones(iZone)%il
        else if (rotDir == 3) then
            ! The zone is one cell wide in the k-direction (first loop in readStructuredCGNS.F90)
            dir1 = zones(izone)%il
            dir2 = zones(izone)%jl
            rotOffset = zones(iZone)%il * zones(iZone)%jl
        end if

        ! Number of nodes on a single plane for this zone
        nNodes = dir1 * dir2

        ! Loop over the nodes on both planes in this zone (mirror and rotated)
        nodeLoop: do iNode = 1, 2 * nNodes
            ii = ii + 1  ! Global index

            if (3 * ii - 2 .lt. istart .or. 3 * ii .gt. iend) then
                ! ii is not in the local ownership range so skip
                cycle nodeLoop
            end if

            ! Determine if the node is on the mirror plane or the rotated plane (1 or 0)
            iPlane = mod((iNode - 1) / rotOffset + 1, 2)
            if (mirrorIdx .eq. 1 .and. iPlane .eq. 1) then
                ! Mirror plane is the first plane and we are on the mirror plane
                jj = ii + rotOffset
            else if (mirrorIdx .eq. 2 .and. iPlane .eq. 0) then
                ! Mirror plane is the second plane and we are on the mirror plane
                jj = ii - rotOffset
            else
                ! This node is on the rot plane so skip
                cycle nodeLoop
            end if

            ! Adjust ii to be in the local indexing
            iiLocal = 3 * ii - istart + 1
            print *, myID, ii, istart, iiLocal

            ! Copy the x, y, z coordinates
            coords(1) = coordsLocal(iiLocal - 2) ! x
            coords(2) = coordsLocal(iiLocal - 1) ! y
            coords(3) = coordsLocal(iiLocal)     ! z

            ! Apply the rotation
            rotCoords = matmul(Mi, coords)

            ! Add the rotated coords to the array of coords to insert
            ! into the global vector
            if (.not. allocated(insertCoords)) then
                allocate (insertCoords(3))
                insertCoords = rotCoords
            else
                insertCoords = [insertCoords, rotCoords]
            end if

            ! Add the index to the array of indices to insert
            ! into the global vector
            if (.not. allocated(insertIdxs)) then
                allocate (insertIdxs(3))
                insertIdxs = (/3 * jj - 2, 3 * jj - 1, 3 * jj/)
            else
                insertIdxs = [insertIdxs, (/3 * jj - 2, 3 * jj - 1, 3 * jj/)]
            end if

        end do nodeLoop

    end do zoneLoop

    ! Insert the rotated coordinates into the global vector
    ! We only do this collective communication once at the end if insertCoords is allocated
    if (allocated(insertCoords)) then
        call VecSetValues(Xv, size(insertCoords), insertIdxs - 1, insertCoords, INSERT_VALUES, ierr)
        call EChk(ierr, __FILE__, __LINE__)
    end if

    ! Need to call VecAssemblyBegin/End
    call VecAssemblyBegin(Xv, ierr)
    call EChk(ierr, __FILE__, __LINE__)
    call VecAssemblyEnd(Xv, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Restore the array to deallocate the coordsLocal pointer
    call VecRestoreArrayF90(Xv, coordsLocal, ierr)
    call EChk(ierr, __FILE__, __LINE__)

end subroutine copyRotateVolumeCoordinates

subroutine copyRotateVolumeCoordinates_d()
    use gridData
    implicit none

    ! Working variables
    integer(kind=intType) :: istart, iend, ierr, iRotStart, ndof, ii, ipatch
    integer(kind=intType) :: a1, b1, a2, idxStart
    real(kind=realType) :: Mi(3, 3), coords_d(3), angle, axis(3)
    real(kind=realType), dimension(:), pointer :: coordsLocal_d, rotCoords_d

    return

    ! do ipatch = 1, size(mirrorPlaneIdxs, 1)
    !     ! Get the necessary axisymm info from the grid data
    !     a1 = mirrorPlaneIdxs(ipatch, 1)  ! Starting index of the mirror plane nodes
    !     b1 = mirrorPlaneIdxs(ipatch, 2) - 1  ! Ending index of the mirror plane nodes
    !     a2 = rotPlaneIdxs(ipatch, 1)  ! Starting index of the rotated nodes
    !     axis = axiSymmAxis  ! Axis of rotation
    !     angle = axiSymmAngle

    !     ! Get the local ownership range
    !     call VecGetOwnershipRange(dXv, istart, iend, ierr)
    !     call EChk(ierr, __FILE__, __LINE__)
    !     istart = istart + 1  ! Convert to 1-based indexing
    !     iend = iend ! Convert to 1-based indexing

    !     call VecGetArrayF90(dXv, coordsLocal_d, ierr)
    !     call EChk(ierr, __FILE__, __LINE__)

    !     if (istart .ge. a1 .and. istart .lt. b1 .and. iend .le. b1) then
    !         ! All local volume nodes lie in [a1, b1]
    !         ndof = iend - istart + 1
    !         idxStart = 1
    !     else if (istart .ge. a1 .and. istart .lt. b1 .and. iend .gt. b1) then
    !         ! istart is in [a1, b1] but iend is not
    !         iend = b1  ! Clip iend to b1
    !         nDof = iend - istart + 1
    !         idxStart = 1
    !     else if (istart .lt. a1 .and. iend .gt. a1 .and. iend .lt. b1) then
    !         ! iend is in [a1, b1] but istart is not
    !         idxStart = a1 - istart + 1
    !         istart = a1  ! Clip istart to a1
    !         nDof = iend - istart + 1
    !     else
    !         ! The local volume nodes are not in [a1, b1]
    !         istart = 0
    !         iend = 0
    !         ndof = 0
    !     end if

    !     ! Get the starting index for the rotated coordinates
    !     iRotStart = a2 + (istart - a1)

    !     call getRotationMatrixAngleAxis(angle, axis, Mi)

    !     allocate (rotCoords_d(ndof))

    !     if (ndof .gt. 0) then
    !         do ii = 1, ndof / 3
    !             ! Extract the x, y, z coordinates
    !             coords_d(1) = coordsLocal_d(idxStart + 3 * ii - 3)
    !             coords_d(2) = coordsLocal_d(idxStart + 3 * ii - 2)
    !             coords_d(3) = coordsLocal_d(idxStart + 3 * ii - 1)

    !             ! Rodrigues rotation formula
    !             rotCoords_d(3 * ii - 2:3 * ii) = matmul(Mi, coords_d)
    !         end do

    !         call VecSetValues(dXv, ndof, (/(ii, ii=iRotStart - 1, iRotStart + nDof - 2, 1)/), &
    !                           rotCoords_d, INSERT_VALUES, ierr)
    !         call EChk(ierr, __FILE__, __LINE__)
    !     end if

    !     ! Need to call VecAssemblyBegin/End
    !     call VecAssemblyBegin(dXv, ierr)
    !     call EChk(ierr, __FILE__, __LINE__)
    !     call VecAssemblyEnd(dXv, ierr)
    !     call EChk(ierr, __FILE__, __LINE__)

    !     call VecRestoreArrayF90(dXv, coordsLocal_d, ierr)
    !     call EChk(ierr, __FILE__, __LINE__)
    ! end do
end subroutine copyRotateVolumeCoordinates_d

subroutine copyRotateVolumeCoordinates_b()
    use gridData
    implicit none

    ! Working variables
    integer(kind=intType) :: istart, iend, ierr, iRotStart, nDofRot, ii, ipatch
    integer(kind=intType) :: a1, b1, a2, b2, iStartLocal, totalDof
    real(kind=realType) :: Mi(3, 3), coords_b(3), angle, axis(3)
    real(kind=realType), dimension(:), pointer :: coordsLocal_b
    real(kind=realType), dimension(:), allocatable :: rotCoords_b

    return

    ! do ipatch = 1, size(mirrorPlaneIdxs, 1)
    !     ! Get the necessary axisymm info from the grid data
    !     a1 = mirrorPlaneIdxs(ipatch, 1)  ! Starting index of the mirror plane nodes
    !     b1 = mirrorPlaneIdxs(ipatch, 2) - 1  ! Ending index of the mirror plane nodes
    !     a2 = rotPlaneIdxs(ipatch, 1)  ! Starting index of the rotated nodes
    !     b2 = rotPlaneIdxs(ipatch, 2) - 1  ! Ending index of the rotated nodes
    !     axis = axiSymmAxis  ! Axis of rotation
    !     angle = axiSymmAngle

    !     ! Get the local ownership range
    !     call VecGetOwnershipRange(dXv, istart, iend, ierr)
    !     call EChk(ierr, __FILE__, __LINE__)
    !     istart = istart + 1  ! Convert to 1-based indexing
    !     iend = iend ! Convert to 1-based indexing

    !     call VecGetArrayF90(dXv, coordsLocal_b, ierr)
    !     call EChk(ierr, __FILE__, __LINE__)

    !     if (istart .ge. a2 .and. iend .le. b2) then
    !         ! All local volume nodes lie in [a2, b2]
    !         nDofRot = iend - istart + 1
    !         iStartLocal = 1
    !     else if (istart .lt. a2 .and. iend .gt. a2 .and. iend .le. b2) then
    !         ! istart is not in [a2, b2] but iend is
    !         iStartLocal = a2 - istart + 1
    !         istart = a2  ! Clip istart to a2
    !         nDofRot = iend - istart + 1
    !     else if (istart .ge. a2 .and. istart .lt. b2 .and. iend .gt. b2) then
    !         ! istart is in [a2, b2] but iend is not
    !         iend = b2  ! Clip iend to b2
    !         nDofRot = iend - istart + 1
    !         iStartLocal = 1
    !     else
    !         istart = 0
    !         iend = 0
    !         nDofRot = 0
    !     end if

    !     ! Get the starting index for the rotated coordinates
    !     iRotStart = a1 + (istart - a2)
    !     allocate (rotCoords_b(nDofRot))

    !     call getRotationMatrixAngleAxis(-angle, axis, Mi)

    !     if (nDofRot .gt. 0) then
    !         do ii = 1, nDofRot / 3
    !             ! Extract the x, y, z coordinates
    !             coords_b(1) = coordsLocal_b(iStartLocal + 3 * ii - 3)
    !             coords_b(2) = coordsLocal_b(iStartLocal + 3 * ii - 2)
    !             coords_b(3) = coordsLocal_b(iStartLocal + 3 * ii - 1)

    !             ! Rodrigues rotation formula
    !             ! Add the rotated seeds to the original seeds
    !             rotCoords_b(3 * ii - 2:3 * ii) = matmul(Mi, coords_b)
    !         end do

    !         ! Set the rotated + original seeds in the correct place
    !         call VecSetValues(dXv, nDofRot, (/(ii, ii=iRotStart - 1, iRotStart + nDofRot - 2)/), &
    !                           rotCoords_b, ADD_VALUES, ierr)
    !         call EChk(ierr, __FILE__, __LINE__)

    !         ! Zero out all seeds on the rotated plane
    !         coordsLocal_b(iStartLocal:iStartLocal + nDofRot - 1) = zero
    !     end if

    !     ! Need to call VecAssemblyBegin/End
    !     call VecAssemblyBegin(dXv, ierr)
    !     call EChk(ierr, __FILE__, __LINE__)
    !     call VecAssemblyEnd(dXv, ierr)
    !     call EChk(ierr, __FILE__, __LINE__)

    !     ! Restore the array to deallocate the coordsLocal_b pointer
    !     call VecRestoreArrayF90(dXv, coordsLocal_b, ierr)
    !     call EChk(ierr, __FILE__, __LINE__)
    ! end do

end subroutine copyRotateVolumeCoordinates_b
