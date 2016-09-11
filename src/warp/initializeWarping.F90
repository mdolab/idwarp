subroutine initializeWarping(pts, uniquePts, link, faceSizes, faceConn, &
     ndofLocal, nUnique, nLink,  nFaceSizes, nFaceConn, restartFile)

  use gridInput
  use communication
  use kd_tree
  use gridData
  implicit none

  ! Input
  real(kind=realType), intent(in) :: pts(ndofLocal), uniquePts(3, nUnique)
  integer(kind=intType), intent(in) :: ndoflocal, nUnique, nLink
  integer(kind=intType), dimension(nFaceSizes), intent(in) :: faceSizes
  integer(kind=intType), dimension(nFaceConn), intent(in) :: faceConn
  integer(kind=intType), dimension(nLink), intent(in) :: link
  integer(kind=intType), intent(in) :: nFaceSizes, nFaceConn
  character*(*) :: restartFile
  ! Working
  integer(kind=intType) :: nNodesTotal, ierr, iStart, iEnd, iProc
  integer(kind=intType) :: i, j, ii, kk
  real(kind=realType), dimension(:), allocatable :: costs, procEndCosts, cumNodesProc
  integer(kind=intType), dimension(:), allocatable :: procSplits, procSplitsLocal
  real(kind=realType), dimension(:), allocatable :: denomenator0Copy
  real(kind=realtype) :: costOffset, totalCost, averageCost, c, tmp, r(3)
  integer(kind=intType) :: nVol, curBin, newBin, newDOFProc, totalDOF
  integer(kind=intType) :: stat
  integer(kind=intType) :: status(MPI_STATUS_SIZE), fileHandle, pointsInFile
  INTEGER(KIND=MPI_OFFSET_KIND) OFFSET
  logical :: recompute, loadFile, saveFile, fileExists, localBadFile, globalBadFile

  ! ----------------------------------------------------------------------
  !   Step 2: Create and set PETSc Vector for the local Nodes and fill
  ! ----------------------------------------------------------------------
  call VecCreateMPI(warp_comm_world, ndofLocal, PETSC_DETERMINE, Xs, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecDuplicate(Xs, dXs, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecGetOwnershipRange(Xs, istart, iend, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  i = 1
  do ii=istart,iend-1
     call VecSetValues(Xs, 1, (/ii/), pts(i), &
          INSERT_VALUES, ierr)
     call EChk(ierr, __FILE__, __LINE__)
     i = i + 1
  end do

  call VecAssemblyBegin(Xs, ierr) 
  call EChk(ierr, __FILE__, __LINE__)
  call VecAssemblyEnd(Xs, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Create a scatter object so everyone can get the full copy Xs
  call VecScatterCreateToAll(Xs, XsToXsLocal, XsLocal, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  Call VecDuplicate(XsLocal, dXsLocal, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  ! For now we just have a single tree....we may have more in the future. 
  allocate(mytrees(1))
  mytrees(1)%tp => myCreateTree(uniquePts, link, faceSizes, faceConn)

  ! We now need to generate the Wi estinate and loadbalance. There are
  ! quite costly and we can save startup time if they are chached by
  ! the user.
  !
  ! There are few different things that can happen here:
  !
  ! 1. No restartFile is supplied so we do the Wi estimate and loadbalance 
  !    as per usual 
  ! 2. Restart file is supplied but doesn't exist, so we save what we found
  ! 3. Load balance file name is supplied and exists. Load and check a few
  !    If they don't match, do regular dry run

  call VecGetArrayF90(commonGridVec, Xv0Ptr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Allocate the denomenator estimate and costs. These may be loaded
  ! from the restart file.
  nVol = size(Xv0Ptr)/3
  allocate(denomenator0(nVol), costs(nVol))
  denomenator0 = zero
  costs = zero

  ! Get the ownership ranges for common grid vector, we need that to
  ! determine where to write things. 
  allocate(cumNodesProc(0:nProc))

  call VecGetOwnershipRanges(commonGridVec, cumNodesProc, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Divide that by 3 since that includes the 3 dof per popint
  cumNodesProc =  cumNodesProc / 3

  if (trim(restartFile) == "") then 
     loadFile = .False.
     recompute = .True.
     saveFile = .False.
  else
     ! See if the file exists
     if (myid == 0) then
        INQUIRE(FILE=trim(restartFile), EXIST=fileExists)
     end if

     call MPI_bcast(fileExists, 1, MPI_LOGICAL, 0, warp_comm_world, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     if (.not. fileExists) then 
        ! Can't load, file doesn't exist. Recompute. 
        loadFile = .False.
        recompute = .True. 
        saveFile = .True.
     else
        ! File exists so we can try to load it. Recompute is still
        ! false, but may be true if the file turns out to be bad.
        loadFile = .True.
        recompute = .False.
        saveFile = .False.
     end if
  end if

  if (loadFile) then 
     if (myid == 0) then 
        print *, 'Loading restart file...'
     end if

     ! Read from the restartFile. We can do this in parallel.
     call mpi_file_open(warp_comm_world, trim(restartFile), &
          MPI_MODE_RDONLY, & 
          MPI_INFO_NULL, fileHandle, ierr) 
     call EChk(ierr, __FILE__, __LINE__)

     if (myid == 0) then 
        ! Write the number of nodes
        call MPI_file_read(fileHandle, pointsInFile, 1, &
             MPI_INTEGER4, MPI_STATUS_IGNORE, ierr)
        call EChk(ierr, __FILE__, __LINE__)
     end if

     call MPI_bcast(pointsInFile, 1, MPI_INTEGER4, 0, warp_comm_world, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     if (pointsInFile /= cumNodesProc(nProc)) then 
        if (myid == 0) then 
           print *, 'Number of points in restart file are different...'
        end if

        ! The file is bad...different number of points. 
        call mpi_file_close(fileHandle, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! Use fortran to delete the file, since mpi_file_delete
        ! doesn't actually work
        if (myid == 0) then 
           open(unit=9, iostat=stat, file=trim(restartFile), status='old')
           if (stat == 0) close(9, status='delete')
        end if

        ! We must now recompute and save
        recompute = .True.
        saveFile = .True.
     else

        ! Continue to read:
        offset = cumNodesProc(myid)*8 + 4
        call mpi_file_seek(fileHandle, offset, MPI_SEEK_SET, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call MPI_file_read(fileHandle, denomenator0, size(denomenator0), &
             MPI_REAL8, status, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        offset = cumNodesProc(nProc)*8 + cumNodesProc(myid)*8 + 4
        call mpi_file_seek(fileHandle, offset, MPI_SEEK_SET, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call MPI_file_read(fileHandle, costs, size(costs), &
             MPI_REAL8, status, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call mpi_file_close(fileHandle, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! It is *still* possible that the data is bad. We will check 
        ! 1/100th of the number of points to see. 

        if (myid == 0) then 
           print *, 'Checking restart file...'
        end if

        localBadFile = .False.

        ! Note that we've swapped the kk loop inside here. This is
        ! such that we can get the full estimate immediately and
        ! verify the correct values
        wiLoop_check: do j=1, nVol, 100
           tmp = zero
           do kk=1, nLoop
              call getMirrorPt(Xv0Ptr(3*j-2:3*j), r, kk)
              call getWiEstimate(mytrees(1)%tp, mytrees(1)%tp%root, r, tmp)
           end do
           ! Do a relative check:
           if (1 - tmp/denomenator0(j) > eps) then 
              localBadFile = .True.
              print *,myid, j, tmp, denomenator0(j), eps
              exit wiLoop_check
           end if
        end do wiLoop_check

        ! Now all reduce to see if anyone thought the file was bad:
        call mpi_AllReduce(localBadFile, globalBadFile, 1, MPI_LOGICAL, &
             MPI_LOR, warp_comm_world, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! It was bad-somehow, so we need to recompute anyway.
        if (globalBadFile) then 
           recompute = .True.
           saveFile = .True.
           if (myid == 0) then 
              print *, 'Restart file not accurate...'
           end if

           ! Also toast the file, since it is useless to us
           if (myid == 0) then 
              open(unit=9, iostat=stat, file=trim(restartFile), status='old')
              if (stat == 0) close(9, status='delete')
           end if
        end if
     end if
  end if

  if (recompute) then 
     if (myid == 0) then 
        print *, 'Computing Denomenator Estimate...'
     end if

     ! Compute the (approx) denomenator. This needs to be done only once.
     do kk=1, nLoop
        wiLoop: do j=1, nVol
           call getMirrorPt(Xv0Ptr(3*j-2:3*j), r, kk)
           call getWiEstimate(mytrees(1)%tp, mytrees(1)%tp%root, r, &
                denomenator0(j))
        end do wiLoop
     end do

     ! For load balacing we have to do a dry run loop that determines
     ! just the number of ops that a mesh warp *would* use. This
     ! integer is scaled by the "brute force" cost, so each individual
     ! cost will be less than 1.0

     if (myid == 0) then 
        print *, 'Load Balancing...'
     end if

     do kk=1, nLoop
        dryRunLoop: do j=1, nVol
           c = zero
           call getMirrorPt(Xv0Ptr(3*j-2:3*j), r, kk)
           call dryRun(mytrees(1)%tp, mytrees(1)%tp%root, r, denomenator0(j), c)
           costs(j) = costs(j) + c / mytrees(1)%tp%n
        end do dryRunLoop
     end do
  end if

  ! Barrier here, otherwise it can look like the saving restart file
  ! takes all the time, in reality, it is waitigng at the collective
  ! write.
  call MPI_barrier(warp_comm_world, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Save the data to the restartFile if required
  if (saveFile) then 
     if (myid == 0) then 
        print *, 'Saving restart file...'
     end if

     ! Write to the file. We can do this in parallel. 
     call mpi_file_open(warp_comm_world, trim(restartFile), &
          MPI_MODE_WRONLY + MPI_MODE_CREATE, & 
          MPI_INFO_NULL, fileHandle, ierr) 
     call EChk(ierr, __FILE__, __LINE__)

     if (myid == 0) then 
        ! Write the number of nodes
        call MPI_file_write(fileHandle, cumNodesProc(nProc), 1, &
             MPI_INTEGER4, MPI_STATUS_IGNORE, ierr)
        call EChk(ierr, __FILE__, __LINE__)
     end if

     offset = cumNodesProc(myid)*8 + 4
     call mpi_file_seek(fileHandle, offset, MPI_SEEK_SET, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     call MPI_file_write(fileHandle, denomenator0, size(denomenator0), &
          MPI_REAL8, MPI_STATUS_IGNORE, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     offset = cumNodesProc(nProc)*8 + cumNodesProc(myid)*8 + 4

     call mpi_file_seek(fileHandle, offset, MPI_SEEK_SET, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     call MPI_file_write(fileHandle, costs, size(costs), &
          MPI_REAL8, MPI_STATUS_IGNORE, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     call mpi_file_close(fileHandle, ierr)
     call EChk(ierr, __FILE__, __LINE__)
  end if

  ! Don't forget to restore arrays!
  call VecRestoreArrayF90(commonGridVec, Xv0Ptr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Now put the costs in cumulative format
  do j=2,nVol
     costs(j) = costs(j) + costs(j-1)
  end do

  ! Now, we gather up the value at the *end* of each of the costs
  ! array and send to everyone
  allocate(procEndCosts(nProc))
  call mpi_allgather(costs(nVol), 1, MPI_REAL8, procEndCosts, 1, MPI_REAL8, &
       warp_comm_world, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Now each proc can "correct" for the cumulative costs by using
  ! the intermediate sums from the other processors
  costOffset = zero
  do iProc=1, myid ! iProc is zero-based but we start on the second proc!
     costOffset = costOffset + procEndCosts(iProc)
  end do

  ! Now just vector sum the costs array
  costs = costs + costOffset

  ! We now have a distributed "costs" array that in cumulative
  ! format has the costs of mesh warping. First thing we will do
  ! make sure all procs know the *precise* total cost. Broadcast
  ! this from the "last" proc
  totalCost = costs(nVol)
  call MPI_bcast(totalCost, 1, MPI_REAL8, nProc-1, warp_comm_world, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! And the average cost:
  averageCost = totalCost / nProc

  ! Next we loop over our own cost array to determine the indicies
  ! in *global* ordering where the "breaks" should be. 
  call vecGetOwnershipRange(commonGridVec, iStart, iEnd, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call vecGetSize(commonGridVec, totalDOF, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  allocate(procSplits(0:nProc), procSplitsLocal(0:nProc))
  procSplitsLocal(:) = 0 ! Do this zero-based ordering!
  curBin = int(costs(1) /  averageCost)
  do i=2, nVol
     newBin = int(costs(i) / averageCost)
     if (newBin > curBin) then
        ! We passed a boundary:
        procSplitsLocal(newBin) = i + iStart/3
        curBin = newBin
     end if
  end do
  ! Make sure the last proc has the right end...
  if (myid == nProc-1) then 
     procSplitsLocal(nProc) = iStart/3 + nVol
  end if
  ! Communicate to all procs 
  call mpi_AllReduce(procSplitsLocal, procSplits, nProc+1, MPI_INTEGER, MPI_MAX, &
       warp_comm_world, ierr)
  call EChk(ierr,__FILE__,__LINE__)
  ! Convert Proc Splits to DOF Format
  procSplits = procSplits * 3
  newDOFProc = procSplits(myid+1) - procSplits(myid)

  ! Now create the final vectors we need:
  call VecCreate(warp_comm_world, Xv, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Set to be be blocked
  call VecSetBlockSize(Xv, 3, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Type and size
  call VecSetType(Xv, "mpi", ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecSetSizes(Xv, newDOFProc, PETSC_DECIDE, ierr)
  call EChk(ierr,__FILE__,__LINE__)
  call VecGetSize(Xv, i, ierr)

  call VecDuplicate(Xv, Xv0, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecDuplicate(Xv, dXv, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Now we know the break-points of the splits, we can create a vec
  ! scatter that converts between the "common" and "warping"
  ! ordering (the repartitioned) ordering. 
  call ISCreateStride(warp_comm_world, newDOFProc, procSplits(myid), 1, IS1, ierr)
  call EChk(ierr, __FILE__, __LINE__)  

  call ISCreateStride (warp_comm_world, newDOFProc, procSplits(myid), 1, IS2, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecScatterCreate(commonGridVec, IS1, Xv, IS2, common_to_warp, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call ISDestroy(IS1, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call ISDestroy(IS2, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! We have to now scatter our nominal grid nodes in
  ! "commonGridVec" to our final ordering. This uses warp_to_common
  ! in REVERSE:
  call VecScatterBegin(common_to_warp, commonGridVec, Xv, &
       INSERT_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecScatterEnd  (common_to_warp, commonGridVec, Xv, &
       INSERT_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! And set Xv0
  call VecCopy(Xv, Xv0, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! We also want to scatter denomenator0 since that was kinda costly
  ! to compute. So we vecPlace it into commonGridVec and take it out
  ! of Xv. However, denomenator is only a third of the size of the
  ! gridVec. So we dump it into a temporary vector, which *is* the
  ! right size and the copy after

  ! We can now also allocate the final space for the denomenator
  allocate(numerator(3, newDOFProc/3))
  allocate(denomenator(newDOFProc/3))
  allocate(denomenator0Copy(nVol*3))
  do i=1,nVol
     denomenator0Copy(3*i-2) = denomenator0(i)
  end do

  ! Clear the current denomenator0 and make the right size.
  deallocate(denomenator0)
  allocate(denomenator0(newDOFProc/3))

  call VecPlaceArray(commonGridVec, denomenator0Copy, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecPlaceArray(Xv, numerator, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Actual scatter
  call VecScatterBegin(common_to_warp, commonGridVec, Xv, &
       INSERT_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecScatterEnd  (common_to_warp, commonGridVec, Xv, &
       INSERT_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! And reset the arrays
  call VecResetArray(commonGridVec, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecResetArray(Xv, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Copy denomenator0 out of 'numerator'
  do i=1,newDOFProc/3
     denomenator0(i) = numerator(1, i)
  end do
  warpMeshDOF = newDOFProc
  commonMeshDOF = nVol*3

  ! Deallocate the memory from this subroutine
  deallocate( cumNodesProc)
  deallocate(costs, procEndCosts, procSplits, procSplitsLocal, denomenator0Copy)

  if (myid == 0) then 
     print *, 'Finished Mesh Initialization.'
  end if
  initializationSet = 1

end subroutine initializeWarping

subroutine setSymmetryPlanes(pts, normals, n)

  use constants
  use gridData

  ! Set the symmetry plane information. 

  real(kind=realType), dimension(3, n), intent(in) :: pts, normals
  integer(kind=intType), intent(in) :: n

  real(kind=realType) :: ptstar(3), r(3)
  integer(kind=intType) :: kk
  ! Loops have be applied recursively so for each symmetry plane, we
  ! multiply the number of iterations by 2. 
  nLoop = 2**(n)

  allocate(symmPts(3, n), symmNormals(3, n))
  symmPts = pts
  symmNormals = normals

end subroutine setSymmetryPlanes

subroutine getMirrorPt(pt, r, kk)

  use gridData

  ! Compute the point reflected point based on the symmetry plane loop
  ! index. We support up to 6 symmetry planes

  ! Input/Output
  real(kind=realType), intent(in) :: pt(3)
  real(kind=realType), intent(out) :: r(3)
  integer(kind=intType), intent(in) ::kk

  ! Working
  real(kind=realType) :: d(3), n(3), dp

  ! Set the input variable to the out
  r = pt

  ! This is the slick implementation. See below for logically expanded
  ! form. 

  do j=1, size(symmPts, 2)
     if (mod(kk-1, 2**j) >= 2**(j-1)) then 
        ! Need to mirror about the jth plane
        d = r - symmPts(:, j)
        n = symmNormals(:, j)
        dp = d(1)*n(1) + d(2)*n(2) + d(3)*n(3)
        r = r - 2*n*dp
     end if
  end do

  ! An explicit implementation would look something like the
  ! following:
  ! if (mod(kk-1, 2) >= 1) then
  !    call mirror(r, symmPts(:, 1), symmNormals(:, 1))
  ! end if

  ! if (mod(kk-1, 4) >= 2) then 
  !    call mirror(r, symmPts(:, 2), symmNormals(:, 2))
  ! end if

  ! if (mod(kk-1, 8) >= 4) then 
  !    call mirror(r, symmPts(:, 3), symmNormals(:, 3))
  ! end if


end subroutine getMirrorPt


subroutine getMirrorNumerator(num, kk)

  ! Mirror the numerator as per necessary based on the loop index
  use gridData

  ! Input/Output
  real(kind=realType), intent(inout) :: num(3)
  integer(kind=intType), intent(in) :: kk

  ! Working
  real(kind=realType) ::  n(3), dp

  do j=1, size(symmPts, 2)
     if (mod(kk-1, 2**j) >= 2**(j-1)) then 
        ! Need to mirror about the jth plane
        n = symmNormals(:, j)
        dp = num(1)*n(1) + num(2)*n(2) + num(3)*n(3)
        num = num - two*dp*n
     end if
  end do

end subroutine getMirrorNumerator

#ifndef USE_COMPLEX
subroutine getMirrorNumerator_b(numb, kk)

  use gridData
  implicit none

  ! Input/Output
  real(kind=realType), intent(inout) :: numb(3)
  integer(kind=intType), intent(in) :: kk

  ! Working
  real(kind=realType) ::  n(3), dp, dpb
  integer(kind=intType) :: branch, j

  DO j=1, size(symmPts, 2)
     IF (MOD(kk - 1, 2**j) .GE. 2**(j-1)) THEN
        CALL PUSHCONTROL1B(1)
     ELSE
        CALL PUSHCONTROL1B(0)
     END IF
  END DO
  DO j=size(symmPts, 2),1,-1
     CALL POPCONTROL1B(branch)
     IF (branch .NE. 0) THEN
        n = symmnormals(:, j)
        dpb = -SUM(n*two*numb)
        numb(1) = numb(1) + n(1)*dpb
        numb(2) = numb(2) + n(2)*dpb
        numb(3) = numb(3) + n(3)*dpb
     END IF
  END DO
end subroutine getMirrorNumerator_b
#endif
subroutine getMirrorNumerator_d(numd, kk)

  use gridData

  implicit none

  ! Input/Output
  real(kind=realType), intent(inout) :: numd(3)
  integer(kind=intType), intent(in) :: kk

  ! Working
  real(kind=realType) ::  n(3), dpd
  integer(kind=intType) :: j

  do j=1, size(symmPts, 2)
     if (mod(kk-1, 2**j) >= 2**(j-1)) then 
        ! Need to mirror about the jth plane
        n = symmnormals(:, j)
        dpd = n(1)*numd(1) + n(2)*numd(2) + n(3)*numd(3)
        numd = numd - two*n*dpd
     END IF
  END DO
  
end subroutine getMirrorNumerator_d
