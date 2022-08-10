subroutine readPlot3d(plot3d_file)
    use gridData
    use gridInput
    use communication
    use CGNSGrid

    implicit none

    ! Input Arguments
    character*(*), intent(in) :: plot3d_file

    ! Working
    integer(kind=intType) :: ierr, iStart
    character(len=12) :: zoneName
    integer(kind=intType) :: i, j, k, ii, iZone, iend, localsize, iProc
    integer(kind=intType) :: nZones, nNodes, dims(3)
    real(kind=realType), dimension(:, :, :), allocatable :: coorX, coorY, coorZ
#ifdef USE_COMPLEX
    complex(kind=realType), dimension(:, :), allocatable :: allNodes, localNodes
#else
    real(kind=realType), dimension(:, :), allocatable :: allNodes, localNodes
#endif
    integer(kind=intType), dimension(:), allocatable :: surfaceNodes, localSurfaceNodes
    integer(kind=intType), dimension(:, :), allocatable :: sizes
    integer(kind=intType) :: status(MPI_STATUS_SIZE)

    if (myid == 0) then
        print *, ' -> Reading PLot3d File: ', plot3d_file

        ! We will be assuming multiblock, unformatted without iblank array.
        open (unit=50, form='unformatted', file=plot3d_file)

        ! Read total number of zones and allocate the zone derived type
        read (50) nZones
        allocate (zones(nZones))

        ! Allocate space for the size array and read.
        allocate (sizes(3, nZones))
        read (50) (sizes(1, i), sizes(2, i), sizes(3, i), i=1, nZones)

        ! Compute the total number of nodes and set the node sizes
        nNodes = 0
        do iZone = 1, nZones
            nNodes = nNodes + sizes(1, iZone) * sizes(2, iZone) * sizes(3, iZone)
            zones(iZone)%il = sizes(1, iZone)
            zones(iZone)%jl = sizes(2, iZone)
            zones(iZone)%kl = sizes(3, iZone)
            write (zoneName, "((a), (I5))") 'Domain.', i
            zones(iZone)%name = zoneName
        end do

        ! Now we know the total number of nodes we can allocate the final
        ! required space and read them in.
        allocate (allNodes(3, nNodes))
        allocate (surfaceNodes(nNodes))
        surfaceNodes = 0

        ii = 0
        zoneLoop: do iZone = 1, nZones
            dims = sizes(:, iZone)
            allocate (coorX(dims(1), dims(2), dims(3)), &
                      coorY(dims(1), dims(2), dims(3)), &
                      coorZ(dims(1), dims(2), dims(3)))

            ! Actual read command
            read (50) &
                (((coorX(i, j, k), i=1, dims(1)), j=1, dims(2)), k=1, dims(3)), &
                (((coorY(i, j, k), i=1, dims(1)), j=1, dims(2)), k=1, dims(3)), &
                (((coorZ(i, j, k), i=1, dims(1)), j=1, dims(2)), k=1, dims(3))

            ! Now interlace the packing
            do k = 1, dims(3)
                do j = 1, dims(2)
                    do i = 1, dims(1)
                        ii = ii + 1
#ifdef USE_COMPLEX
                        allNodes(1, ii) = cmplx(coorX(i, j, k), 0.0)
                        allNodes(2, ii) = cmplx(coorY(i, j, k), 0.0)
                        allNodes(3, ii) = cmplx(coorZ(i, j, k), 0.0)
#else
                        allNodes(1, ii) = coorX(i, j, k)
                        allNodes(2, ii) = coorY(i, j, k)
                        allNodes(3, ii) = coorZ(i, j, k)
#endif
                    end do
                end do
            end do

            ! And free the temporary arrays
            deallocate (coorX, coorY, coorZ)
        end do zoneLoop

        ! Close the file
        close (50)
    end if

    ! Communicate the total number of nodes to everyone
    call MPI_bcast(nNodes, 1, MPI_INTEGER, 0, warp_comm_world, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    rootProc2: if (myid == 0) then

        istart = 0 + 1
        iend = int(nNodes * (one / nProc))
        localSize = iend - istart + 1
        allocate (localNodes(3, localSize))
        allocate (localSurfaceNodes(localSize))

        ! Just copy for the root proc:
        localNodes(:, :) = allNodes(:, 1:localSize)
        localSurfaceNodes(:) = surfaceNodes(1:localSize)

        ! Loop over the remainer of the procs and send
        do iProc = 1, nProc - 1
            istart = int(nNodes * (dble(iProc) / nProc)) + 1
            iend = int(nNodes * (dble(iProc + 1) / nProc))
            localSize = iend - istart + 1
#ifdef USE_COMPLEX
            call MPI_Send(allNodes(:, iStart), localSize * 3, MPI_COMPLEX16, iProc, &
                          11, warp_comm_world, ierr)
#else
            call MPI_Send(allNodes(:, iStart), localSize * 3, MPI_REAL8, iProc, &
                          11, warp_comm_world, ierr)
#endif
            call EChk(ierr, __FILE__, __LINE__)

            call MPI_Send(surfaceNodes(iStart), localSize, MPI_INTEGER4, iProc, &
                          12, warp_comm_world, ierr)
            call EChk(ierr, __FILE__, __LINE__)
        end do
        deallocate (allNodes)
        deallocate (surfaceNodes)
    else

        istart = int(nNodes * (dble(myid) / nProc)) + 1
        iend = int(nNodes * (dble(myid + 1) / nProc))
        localSize = iend - istart + 1

        allocate (localNodes(3, localSize))
        allocate (localSurfaceNodes(localSize))

        ! Receive on all the other procs:
#ifdef USE_COMPLEX
        call MPI_recv(localNodes, 3 * localSize, MPI_COMPLEX16, 0, 11, &
                      warp_comm_world, status, ierr)
#else
        call MPI_recv(localNodes, 3 * localSize, MPI_REAL8, 0, 11, &
                      warp_comm_world, status, ierr)
#endif
        call EChk(ierr, __FILE__, __LINE__)
        call MPI_recv(localSurfaceNodes, localSize, MPI_INTEGER4, 0, 12, &
                      warp_comm_world, status, ierr)
        call EChk(ierr, __FILE__, __LINE__)
    end if rootProc2

    ! All we are doing to do here is to create the Xv vector --- which
    ! is done via a call to createCommonGrid

    call createCommonGrid(localNodes, size(localNodes, 2))
    deallocate (localNodes, localSurfaceNodes)

end subroutine readPlot3d

subroutine readPlot3dSurface(plot3d_file)

    ! Read the plot3d surface file that will be used to define
    use communication
    use plot3dSurface

    implicit none

    ! Input Arguments
    character*(*), intent(in) :: plot3d_file

    ! Working
    integer(kind=intType) :: i, j, k, ii, jj, iZone, dims(2)
    integer(kind=intType) :: nZones, nNodes, nElems
    real(kind=realType), dimension(:, :), allocatable :: coorX, coorY, coorZ
    integer(kind=intType), dimension(:, :), allocatable :: sizes

    if (myid == 0) then
        print *, ' -> Reading Plot3d Surface File: ', plot3d_file

        ! We will be assuming multiblock, unformatted without iblank array.
        open (unit=50, form='unformatted', file=plot3d_file)

        ! Read total number of zones and allocate the zone derived type
        read (50) nZones

        ! Allocate space for the size array and read.
        allocate (sizes(2, nZones))
        read (50) (sizes(1, i), sizes(2, i), i=1, nZones)

        ! Compute the total number of nodes and set the node sizes
        nNodes = 0
        nElems = 0
        do iZone = 1, nZones
            nNodes = nNodes + sizes(1, iZone) * sizes(2, iZone)
            nElems = nElems + (sizes(1, iZone) - 1) * (sizes(2, iZone) - 1)
        end do

        ! Now we know the total number of nodes we can allocate the final
        ! required space and read them in.
        allocate (pts(3, nNodes))
        allocate (conn(4, nElems))

        ii = 0
        jj = 0

        zoneLoop: do iZone = 1, nZones
            dims = sizes(:, iZone)
            allocate (coorX(dims(1), dims(2)), &
                      coorY(dims(1), dims(2)), &
                      coorZ(dims(1), dims(2)))

            ! Actual node read command
            read (50) &
                ((coorX(i, j), i=1, dims(1)), j=1, dims(2)), &
                ((coorY(i, j), i=1, dims(1)), j=1, dims(2)), &
                ((coorZ(i, j), i=1, dims(1)), j=1, dims(2))

            ! Now generate the connectivity. We are going to have to
            ! assume here that the orientation of the surface is correct
            ! since we have no way of checking it.  Do this before the
            ! actual interlacing that we ii is un-incremented and it
            ! always has the correct offset.

            do j = 0, dims(2) - 2
                do i = 0, dims(1) - 2
                    jj = jj + 1
                    conn(1, jj) = ii + (j) * dims(1) + i + 1     ! n1
                    conn(2, jj) = ii + (j) * dims(1) + i + 1 + 1 ! n2
                    conn(3, jj) = ii + (j + 1) * dims(1) + i + 1 + 1 ! n3
                    conn(4, jj) = ii + (j + 1) * dims(1) + i + 1     ! n4
                end do
            end do

            ! Now interlace the packing
            do j = 1, dims(2)
                do i = 1, dims(1)
                    ii = ii + 1
#ifdef USE_COMPLEX
                    pts(1, ii) = cmplx(coorX(i, j), 0.0)
                    pts(2, ii) = cmplx(coorY(i, j), 0.0)
                    pts(3, ii) = cmplx(coorZ(i, j), 0.0)
#else
                    pts(1, ii) = coorX(i, j)
                    pts(2, ii) = coorY(i, j)
                    pts(3, ii) = coorZ(i, j)
#endif
                end do
            end do

            ! And free the temporary arrays
            deallocate (coorX, coorY, coorZ)
        end do zoneLoop

        ! Close the file
        close (50)
        deallocate (sizes)
    end if
end subroutine readPlot3dSurface
