subroutine writePlot3d(plot3d_file)

    use precision
    use communication
    use CGNSGrid
    use gridData
    implicit none

    ! Input Arguments
    character*(*), intent(in) :: plot3d_file

    ! Working variables
    integer(kind=intType) :: i, j, jj, offset, nNodes, iZone, ierr
    real(kind=realType), dimension(:), allocatable :: coorX, coorY, coorZ
    real(kind=realType), pointer, dimension(:) :: xx
    integer(kind=intType), dimension(:, :), allocatable :: sizes

    ! This isn't technically scalable...we will dump the entire grid
    ! onto the root proc and write there.
    call VecScatterCreateToZero(Xv, XvToLocal, Xvlocal, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call VecScatterBegin(XvToLocal, Xv, XvLocal, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call VecScatterEnd(XvToLocal, Xv, XvLocal, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Only do writing on root proc:
    rootProc: if (myid == 0) then

        call VecGetArrayF90(XvLocal, xx, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! Open the new file
        open (unit=50, form='unformatted', file=plot3d_file)

        ! Write the plot3d header
        write (50) size(zones)

        ! Gather the zone sizes
        allocate (sizes(3, size(zones)))
        do iZone = 1, size(zones)
            sizes(1, iZone) = zones(iZone)%il
            sizes(2, iZone) = zones(iZone)%jl
            sizes(3, iZone) = zones(iZone)%kl
        end do

        ! Write the zone sizes
        write (50) (sizes(1, i), sizes(2, i), sizes(3, i), i=1, size(zones))

        offset = 0
        do iZone = 1, size(zones)

            nNodes = zones(iZone)%il * zones(iZone)%jl * zones(iZone)%kl

            ! Allocate arrays for X,Y and Z.Note that we CANNOT just use
            ! matrix notation and stride directly from xx since that will
            ! cause a stack copy to be creater which will make ifort barf,
            ! since they put temps are the stack by default *shakes fist*
            allocate (coorX(nNodes), coorY(nNodes), coorZ(nNodes))

            ! Extract the coordinates
            do j = 1, nNodes
                jj = offset + j
                coorX(j) = xx(3 * jj - 2)
                coorY(j) = xx(3 * jj - 1)
                coorZ(j) = xx(3 * jj)
            end do

            offset = offset + nNodes

            ! Actual write command
            write (50) &
                (coorX(i), i=1, nNodes), &
                (coorY(i), i=1, nNodes), &
                (coorZ(i), i=1, nNodes)

            deallocate (coorX, coorY, coorZ)
        end do

        call VecRestoreArrayF90(XvLocal, xx, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! Close the output file
        close (50)

        ! Free the sizes array only on root proc
        deallocate (sizes)
    end if rootProc

    call vecDestroy(XvLocal, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call VecScatterDestroy(XvToLocal, ierr)
    call EChk(ierr, __FILE__, __LINE__)

end subroutine writePlot3d

