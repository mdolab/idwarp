subroutine readStructuredCGNS(cg)
    use gridData
    use gridInput
    use communication
    use CGNSGrid

    implicit none

    ! Input Arguments
    integer(kind=intType) :: cg

    ! Working
    integer(kind=intType) :: i, j, k, ii, istart, iend, localsize, iProc, iZone, offset
    integer(kind=intType) :: ierr, base, nZones, nNodes
    integer(kind=cgsize_t) :: dims(9), npts, tmp(3, 2)
    integer(kind=intType) :: start(3)
    character(len=maxCGNSNameLen) :: bocoName, famName, coorName(3), zoneName
    integer(kind=intType) :: nbocos, bocotype, boco
    integer(kind=intType) :: ptset_type, normalIndex(3), datatype, ndataset, pts(3, 2)
    integer(kind=intType) :: coorDataType(3)
    integer(kind=cgsize_t) :: normalListFlag
    real(kind=8) :: data_double(6)
    real(kind=realType), dimension(:, :, :), allocatable :: coorX, coorY, coorZ
#ifdef USE_COMPLEX
    complex(kind=realType), dimension(:, :), allocatable :: allNodes, localNodes
#else
    real(kind=realType), dimension(:, :), allocatable :: allNodes, localNodes
#endif
    integer(kind=intType), dimension(:), allocatable :: surfaceNodes, localSurfaceNodes
    integer(kind=intType), dimension(:, :), allocatable :: sizes
    integer(kind=intType) :: nSurf, nodeCount, nConn, ni, nj, nx, ny, nz, bocoIdx
    integer(kind=intType) :: status(MPI_STATUS_SIZE)
    integer(kind=intType) :: nTotalBocos
    logical :: lowFace

    ! Only do reading on root proc:
    rootProc: if (myid == 0) then

        ! Always the first base
        base = 1_intType

        call cg_nzones_f(cg, base, nZones, ierr); 
        if (ierr .eq. CG_ERROR) call cg_error_exit_f
        allocate (zones(nZones))

        ! Now count up the total number of nodes
        nNodes = 0
        allocate (sizes(3, nZones))
        do i = 1, nZones
            call cg_zone_read_f(cg, base, i, zoneName, dims, ierr)
            sizes(:, i) = dims(1:3)
            nNodes = nNodes + dims(1) * dims(2) * dims(3)

            ! Nullify section pointers since we won't have any for structured mesh
            nullify (zones(i)%sections)
        end do

        ! Now we know the total number of nodes we can allocate the final
        ! required space and read them in.
        allocate (allNodes(3, nNodes))
        allocate (surfaceNodes(nNodes))
        surfaceNodes = 0

        ! We make two loops over the zones....the first *just* reads the
        ! BC info and stores it. The second, actually reads the
        ! nodes.

        nSurf = 0
        nConn = 0
        nTotalBocos = 0
        bocoIdx = 0
        zoneLoop1: do iZone = 1, nZones
            call cg_zone_read_f(cg, base, iZone, zoneName, dims, ierr)
            if (ierr .eq. CG_ERROR) call cg_error_exit_f
            zones(iZone)%il = dims(1)
            zones(iZone)%jl = dims(2)
            zones(iZone)%kl = dims(3)
            zones(iZone)%name = zoneName

            call cg_nbocos_f(cg, base, iZone, nbocos, ierr)
            if (ierr .eq. ERROR) call cg_error_exit_f
            allocate (zones(iZone)%bocos(nbocos))
            nTotalBocos = nTotalBocos + nbocos

            bocoLoop1: do boco = 1, nbocos

                ! Nullify elemPtr, elemConn and elemNodes since we won't have any
                nullify (zones(iZone)%bocos(boco)%elemPtr, &
                         zones(iZone)%bocos(boco)%elemConn, &
                         zones(iZone)%bocos(boco)%elemNodes)

                call cg_boco_info_f(cg, base, iZone, boco, boconame, bocotype, &
                                    ptset_type, npts, NormalIndex, NormalListFlag, datatype, &
                                    ndataset, ierr)
                if (ierr .eq. CG_ERROR) call cg_error_exit_f

                call cg_boco_read_f(cg, base, iZone, boco, tmp, data_double, ierr)
                if (ierr .eq. CG_ERROR) call cg_error_exit_f
                pts = int(tmp, intType)
                ! Save the boco info
                zones(iZone)%bocos(boco)%name = boconame
                zones(iZone)%bocos(boco)%type = bocoType
                zones(iZone)%bocos(boco)%ptRange = pts
                zones(iZone)%bocos(boco)%family = ""

                ! Read family  name if possible
                call cg_goto_f(cg, base, ierr, "Zone_t", iZone, "ZoneBC_t", 1, "BC_t", boco, "end")
                if (ierr == 0) then ! Node exists
                    call cg_famname_read_f(famName, ierr)

                    ! Famname exists only if there was no error:
                    if (ierr .ne. CG_ERROR) then
                        zones(iZone)%bocos(boco)%family = famName
                    else
                        ! Set a default name based on the type. We have a
                        ! map table setup for that
                        zones(iZone)%bocos(boco)%family = defaultFamName(bocoType)
                    end if

                    ! Actual number of nodes on patch
                    nx = pts(1, 2) - pts(1, 1) + 1
                    ny = pts(2, 2) - pts(2, 1) + 1
                    nz = pts(3, 2) - pts(3, 1) + 1

                    if (pts(1, 1) == pts(1, 2)) then
                        nSurf = nSurf + ny * nz
                        nConn = nConn + (ny - 1) * (nz - 1)
                    else if (pts(2, 1) == pts(2, 2)) then
                        nSurf = nSurf + nx * nz
                        nConn = nConn + (nx - 1) * (nz - 1)
                    else if (pts(3, 1) == pts(3, 2)) then
                        nSurf = nSurf + nx * ny
                        nConn = nConn + (nx - 1) * (ny - 1)
                    end if
                end if
            end do bocoLoop1
        end do zoneLoop1

        ! Allocate space for storage of all the surface nodes and the
        ! connectivity
        allocate (surfacePoints(3 * nSurf))
        allocate (surfaceConn(4 * nConn))

        ! Allocate the flattened BC types storage
        allocate (bocoTypes(nTotalBocos))

        ! Reset the counters here
        nSurf = 0
        offSet = 0
        nConn = 0
        nodeCount = 0
        ii = 0
        zoneLoop2: do iZone = 1, nZones
            dims(1) = zones(iZone)%il
            dims(2) = zones(iZone)%jl
            dims(3) = zones(iZone)%kl

            start = (/1, 1, 1/)
            allocate (coorX(dims(1), dims(2), dims(3)), &
                      coorY(dims(1), dims(2), dims(3)), &
                      coorZ(dims(1), dims(2), dims(3)))

            ! Check we have doubles in the gridFile
            do i = 1, 3
                call cg_coord_info_f(cg, base, iZone, i, coorDataType(i), coorName(i), ierr)
                if (ierr .eq. CG_ERROR) call cg_error_exit_f
                if (coorDataType(i) /= RealDouble) then
                    write (*, "((a), (I4), (a), (I4), (a), (I1), (a))") 'Error: Coordinates for base ', &
                        base, ', zone ', iZone, ', coordinate ', i, &
                        ', are are not double precision. All coordiante arrays must be double precision.'
                    stop
                end if
            end do

          call cg_coord_read_f(cg, base, iZone, "CoordinateX", RealDouble, int(start, cgsize_t), dims(1:3), coorX, ierr)
            if (ierr .eq. CG_ERROR) call cg_error_exit_f
          call cg_coord_read_f(cg, base, iZone, "CoordinateY", RealDouble, int(start, cgsize_t), dims(1:3), coorY, ierr)
            if (ierr .eq. CG_ERROR) call cg_error_exit_f
          call cg_coord_read_f(cg, base, iZone, "CoordinateZ", RealDouble, int(start, cgsize_t), dims(1:3), coorZ, ierr)
            if (ierr .eq. CG_ERROR) call cg_error_exit_f

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

            bocoLoop2: do boco = 1, size(zones(iZone)%bocos)
                bocoIdx = bocoIdx + 1
                bocoType = zones(iZone)%bocos(boco)%type
                bocoTypes(bocoIdx) = bocoType
                pts = zones(iZone)%bocos(boco)%ptRange

                ! Flag the surface nodes
                do k = pts(3, 1), pts(3, 2)
                    do j = pts(2, 1), pts(2, 2)
                        do i = pts(1, 1), pts(1, 2)
                            ! Surface nodes in global ordering
                            surfaceNodes(offset + (k - 1) * dims(1) * dims(2) + &
                                         (j - 1) * dims(1) + &
                                         i) = 1

                            ! List of all surface nodes
                            surfacePoints(3 * nSurf + 1) = coorX(i, j, k)
                            surfacePoints(3 * nSurf + 2) = coorY(i, j, k)
                            surfacePoints(3 * nSurf + 3) = coorZ(i, j, k)
                            nSurf = nSurf + 1
                        end do
                    end do
                end do

                lowFace = .False.
                ! Determine generic face size
                if (pts(1, 1) == pts(1, 2)) then ! iMin/iMax
                    ni = pts(2, 2) - pts(2, 1) + 1
                    nj = pts(3, 2) - pts(3, 1) + 1
                    if (pts(1, 1) == 1) then
                        lowFace = .True.
                    end if

                else if (pts(2, 1) == pts(2, 2)) then !jMin/jMax
                    ni = pts(1, 2) - pts(1, 1) + 1
                    nj = pts(3, 2) - pts(3, 1) + 1
                    if (pts(2, 1) == 1) then
                        lowFace = .True.
                    end if

                else ! kMin/kMax
                    ni = pts(1, 2) - pts(1, 1) + 1
                    nj = pts(2, 2) - pts(2, 1) + 1

                    if (pts(3, 1) == 1) then
                        lowFace = .True.
                    end if
                end if

                ! Loop over generic face size...We are doing 1 based
                ! ordering. If it is low face normal ordering:
                !
                ! i, j+1 +-----+ i+1, j+1
                !   n4   |     | n3
                !        +-----+
                !       i,j    i+1, j
                !       n1     n2
                !
                ! Otherwise, we flip the ordering to be n1, n4, n3, n2
                if (lowFace) then
                    do j = 0, nj - 2
                        do i = 0, ni - 2
                            surfaceConn(4 * nConn + 1) = nodeCount + (j) * ni + i + 1     ! n1
                            surfaceConn(4 * nConn + 2) = nodeCount + (j) * ni + i + 1 + 1 ! n2
                            surfaceConn(4 * nConn + 3) = nodeCount + (j + 1) * ni + i + 1 + 1 ! n3
                            surfaceConn(4 * nConn + 4) = nodeCount + (j + 1) * ni + i + 1     ! n4
                            nConn = nConn + 1
                        end do
                    end do
                else ! Flip the orientation:
                    do j = 0, nj - 2
                        do i = 0, ni - 2
                            surfaceConn(4 * nConn + 1) = nodeCount + (j) * ni + i + 1     ! n1
                            surfaceConn(4 * nConn + 2) = nodeCount + (j + 1) * ni + i + 1     ! n4
                            surfaceConn(4 * nConn + 3) = nodeCount + (j + 1) * ni + i + 1 + 1 ! n3
                            surfaceConn(4 * nConn + 4) = nodeCount + (j) * ni + i + 1 + 1 ! n2
                            nConn = nConn + 1
                        end do
                    end do
                end if
                nodeCount = nodeCount + ni * nj

            end do bocoLoop2

            offset = offset + dims(1) * dims(2) * dims(3)

            ! And free the temporary arrays
            deallocate (coorX, coorY, coorZ)

        end do zoneLoop2

        ! All finished with the CGNS file.
        call cg_close_f(cg, ierr)
        if (ierr .eq. CG_ERROR) call cg_error_exit_f
        deallocate (sizes)
    else
        !Allocate these to zero so we can just blindly dealloc later
        allocate (surfacePoints(0), surfaceConn(0))
    end if rootProc

    ! Communicate the total number of nodes to everyone
    call MPI_bcast(nNodes, 1, MPI_INTEGER, 0, warp_comm_world, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! For axisymm cases we need to communicate some of the zone data to every proc
    axisymmetric: if (axisymm) then
        ! Total number of zones
        call MPI_bcast(nZones, 1, MPI_INTEGER, 0, warp_comm_world, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! Allocate the zones
        if (myid /= 0) then
            allocate (zones(nZones))
        end if

        ! Communicate the zone data
        do i = 1, nZones
            ! Communicate the zone data
            call MPI_bcast(zones(i)%il, 1, MPI_INTEGER, 0, warp_comm_world, ierr)
            call EChk(ierr, __FILE__, __LINE__)
            call MPI_bcast(zones(i)%jl, 1, MPI_INTEGER, 0, warp_comm_world, ierr)
            call EChk(ierr, __FILE__, __LINE__)
            call MPI_bcast(zones(i)%kl, 1, MPI_INTEGER, 0, warp_comm_world, ierr)
            call EChk(ierr, __FILE__, __LINE__)
            call MPI_bcast(zones(i)%name, maxCGNSNameLen, MPI_CHARACTER, 0, warp_comm_world, ierr)
            call EChk(ierr, __FILE__, __LINE__)
            call MPI_bcast(zones(i)%nVertices, 1, MPI_INTEGER, 0, warp_comm_world, ierr)
            call EChk(ierr, __FILE__, __LINE__)
            call MPI_bcast(zones(i)%nElements, 1, MPI_INTEGER, 0, warp_comm_world, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            ! Get the number of bocos for this zone from the root proc
            if (myID == 0) nbocos = size(zones(i)%bocos)

            ! Communicate the number of bocos
            call MPI_bcast(nbocos, 1, MPI_INTEGER, 0, warp_comm_world, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            if (myid /= 0) then
                allocate (zones(i)%bocos(nbocos))

                ! Nullify section pointers since we won't have any for structured mesh
                nullify (zones(i)%sections)
            end if

            do j = 1, nbocos
                ! Communicate all boco data
                call MPI_bcast(zones(i)%bocos(j)%name, maxCGNSNameLen, MPI_CHARACTER, 0, warp_comm_world, ierr)
                call EChk(ierr, __FILE__, __LINE__)
                call MPI_bcast(zones(i)%bocos(j)%type, 1, MPI_INTEGER, 0, warp_comm_world, ierr)
                call EChk(ierr, __FILE__, __LINE__)
                call MPI_bcast(zones(i)%bocos(j)%ptRange, 6, MPI_INTEGER, 0, warp_comm_world, ierr)
                call EChk(ierr, __FILE__, __LINE__)
                call MPI_bcast(zones(i)%bocos(j)%family, maxCGNSNameLen, MPI_CHARACTER, 0, warp_comm_world, ierr)
                call EChk(ierr, __FILE__, __LINE__)
                call MPI_bcast(zones(i)%bocos(j)%nBCElem, 1, MPI_INTEGER, 0, warp_comm_world, ierr)
                call EChk(ierr, __FILE__, __LINE__)
                call MPI_bcast(zones(i)%bocos(j)%nBCNodes, 1, MPI_INTEGER, 0, warp_comm_world, ierr)
                call EChk(ierr, __FILE__, __LINE__)

                ! Allocate empty data so we don't get errors trying to deallocate later
                if (myid /= 0) then
                    allocate (zones(i)%bocos(j)%BCElements(0))
                    nullify (zones(i)%bocos(j)%elemPtr)
                    nullify (zones(i)%bocos(j)%elemConn)
                    allocate (zones(i)%bocos(j)%elemNodes(0, 0))
                end if
            end do
        end do
    end if axisymmetric

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

end subroutine readStructuredCGNS

subroutine checkInFamilyList(familyList, famName, nFamily, index)
    use precision
    use constants
    implicit none
    character(len=32), dimension(maxFamilies) :: familyList
    character(len=32) :: famName
    integer(kind=intType) :: nFamily, i, index

    index = 0_intType
    famLoop: do i = 1, nFamily
        if (famName == familyList(i)) then
            index = i
            exit famLoop
        end if
    end do famLoop

end subroutine checkInFamilyList
