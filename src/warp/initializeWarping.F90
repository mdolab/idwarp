subroutine initializeWarping(pts, ndof, faceSizesLocal, faceConnLocal, nFaceSizesLocal, &
     nFaceConnLocal)
  
  use gridInput
  use communication
  use kd_tree
  use gridData
  implicit none
 
  ! Input
  real(kind=realType), intent(in) :: pts(ndof)
  integer(kind=intType), intent(in) :: ndof
  integer(kind=intType), dimension(nFaceSizesLocal), intent(in) :: faceSizesLocal
  integer(kind=intTYpe), dimension(nFaceConnLocal), intent(in) :: faceConnLocal
  integer(kind=intType), intent(in) :: nFaceSizesLocal, nFaceConnLocal

  ! Working
  integer(kind=intType) :: nNodesTotal, ierr, iStart, iEnd, iProc
  integer(kind=intType) :: i, j, k, ii, kk, n, ind
  integer(kind=intType), dimension(:), allocatable :: nNodesProc, cumNodesProc
  real(kind=realType),dimension(:,:), allocatable :: allNodes
  real(kind=realType), dimension(:), allocatable :: costs, procEndCosts
  integer(kind=intType), dimension(:), allocatable :: procSplits, procSplitsLocal
  integer(kind=intType), dimension(:), allocatable :: tempInt
  real(kind=realType), dimension(:), allocatable :: denomenator0Copy
  integer(kind=intType), pointer :: indices(:)
  integer(kind=intType), dimension(:), allocatable :: dXsIndices
  real(kind=realtype) :: costOffset, totalCost, averageCost, c, fact(3, 2)
  integer(kind=intType) :: nVol, curBin, newBin, newDOFProc, totalDOF, nFace, nLoop
  real(Kind=realType) :: dists(1), pt(3)
  integer(kind=intType) :: resInd(1), surfID

  ! ----------------------------------------------------------------------
  !   Step 1: Communicate all points with every other proc
  ! ----------------------------------------------------------------------
  allocate(nNodesProc(nProc), cumNodesProc(0:nProc))
  nNodesProc(:) = 0_intType
  
  call mpi_allgather(ndof/3, 1, MPI_INTEGER, nNodesProc, 1, mpi_integer4, &
       warp_comm_world, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Sum and Allocate receive displ offsets
  cumNodesProc(0) = 0_intType
  nNodesTotal = 0
  do i=1, nProc
     nNodesTotal = nNodesTotal + nNodesProc(i)
     cumNodesProc(i) = cumNodesProc(i-1) + nNodesProc(i)
  end do
  nNodesTotal = nNodesTotal
  allocate(allNodes(3, nNodesTotal))
  allNodes = zero ! Make valgrind happy

  call mpi_allgatherv(&
       pts, ndof, MPI_REAL8, &
       allNodes, nNodesProc*3, cumNodesProc*3, MPI_REAL8, &
       warp_comm_world, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! ----------------------------------------------------------------------
  !   Step 2: Create and set PETSc Vector for the local Nodes and fill
  ! ----------------------------------------------------------------------
  call VecCreateMPI(warp_comm_world, ndof, PETSC_DETERMINE, Xs, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecDuplicate(Xs, dXs, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecGetOwnershipRange(Xs, istart, iend, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  i = 1
  do ii=istart,iend-1
#ifdef USE_COMPLEX
     call VecSetValues(Xs, 1, (/ii/), cmplx(pts(i), 0.0), &
          INSERT_VALUES, ierr)
#else
     call VecSetValues(Xs, 1, (/ii/), pts(i), &
          INSERT_VALUES, ierr)
#endif
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

  ! For now we just have a single tree....we may have more in the future. 
  allocate(mytrees(1))
  mytrees(1)%tp => myCreateTree(allNodes, cumNodesProc, faceSizesLocal, faceConnLocal)
  
  ! Compute the approximate denomenator:
  call VecGetArrayF90(commonGridVec, Xv0Ptr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  if (myid == 0) then 
     print *, 'Computing Denomenator Estimate...'
  end if

  nVol = size(Xv0Ptr)/3
  allocate(denomenator0(nVol))

  ! Compute the (approx) denomenator. This needs to be done only once.
  denomenator0 = zero
  call getLoopFact(nLoop, fact)
  do kk=1, nLoop
     wiLoop: do j=1, nVol
        call getWiEstimate(mytrees(1)%tp, mytrees(1)%tp%root, Xv0Ptr(3*j-2:3*j)*fact(:, kk), &
        denomenator0(j))
     end do wiLoop
  end do

  ! If we want to load balance, which generally should be done, we
  ! have to do a dry run loop that determines just the number of
  ! ops that a mesh warp *would* use. This integer is scaled by the
  ! "brute force" cost, so each individual cost will be less than
  ! 1.0
  if (myid == 0) then 
     print *, 'Load Balancing...'
  end if
  allocate(costs(nVol))
  costs = zero
  do kk=1, nLoop
     dryRunLoop: do j=1, nVol
        c = zero
        call dryRun(mytrees(1)%tp, mytrees(1)%tp%root, Xv0Ptr(3*j-2:3*j)*fact(:, kk), &
             denomenator0(j), c)
        costs(j) = costs(j) + c / mytrees(1)%tp%n
     end do dryRunLoop
  end do
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
  deallocate(nNodesProc, cumNodesProc, allNodes)
  deallocate(costs, procEndCosts, procSplits, procSplitsLocal, denomenator0Copy)

  ! --------------- We may need this at some point.... -------------
  ! if (myid == 0) then 
  !    print *, 'Computing wall distance ...'
  ! end if
  ! allocate(d2wall(newDOFProc/3))
  ! call VecGetArrayF90(Xv0, Xv0Ptr, ierr)
  ! call EChk(ierr,__FILE__,__LINE__)

  ! do j=1,warpMeshDOF/3
  !    call n_nearest_to(mytrees(1)%tp, Xv0Ptr(3*j-2:3*j), 1, resInd, dists)
  !    d2wall(j) = dists(1)
  ! end do
  ! call VecRestoreArrayF90(Xv0, Xv0Ptr, ierr)
  ! call EChk(ierr,__FILE__,__LINE__)
  ! -------------------------------------------------------------------

  ! ------------ This is optional -----------

  ! We could perform a test warp to get the "exact" denomenator that
  ! will result from an actual warp. This most likely isn't necessary
  ! however.

  ! if (myid == 0) then
  !    print *, 'Performing Test Warp ...'
  ! end if
  ! call warpMesh()

  ! ! Copy in the real denomenator
  ! denomenator0 = denomenator
  ! -------------------------------------------

  ! One last thing...we need to compute a mapping from the volume
  ! nodes that happen to be on the surface to the actual surface
  ! nodes. This is required for the fake mesh warp derivatives. We we
  ! don't have any addiitonal information we will use the tree do to a
  ! spatial search. This should be reasonably fast since we only are
  ! checking volume nodes we *know* are already on the surface. This
  ! is similar to what is done in pyWarp. 

  if (myid == 0) then 
     print *, 'Performing surface matching ...'
  end if
  allocate(dXsIndices(size(wallIndices)))

  ! We happen to know that all of the indices in wallNodesInXv are
  ! local...so we can use VecGetValues()
  call VecGetOwnershipRange(CommonGridVec, istart, iend, ierr)
  call EChk(ierr, __FILE__, __LINE__)


  do i=1, size(wallIndices)/3
     call VecGetValues(CommonGridVec, 3, (/wallIndices(3*i-2), wallIndices(3*i-1), wallIndices(3*i)/), &
          pt, ierr)
     call EChk(ierr, __FILE__, __LINE__)
    ! print *,myid, i, pt
     ! Now search for that node:
     surfID = pt_in_tree(mytrees(1)%tp, pt)
     if (surfID /= 0) then 
        dXsIndices(3*i-2) = (surfID-1)*3
        dXsIndices(3*i-1) = (surfID-1)*3 + 1
        dXsIndices(3*i  ) = (surfID-1)*3 + 2
     else
        print *,'Error: Could not find a surface node. Are you actually using the same mesh??'
        print *, 'The Bad point is:', pt
        stop
     end if
  end do
  
  call ISCreateGeneral(warp_comm_world, size(wallIndices), wallIndices, &
       PETSC_COPY_VALUES, IS1, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Now create the scatter:
  call ISCreateGeneral(warp_comm_world, size(dXsIndices), dXsIndices, &
       PETSC_COPY_VALUES, IS2, ierr)

  call VecScatterCreate(commonGridVec, IS1, dXs, IS2, common_to_dxs, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call ISDestroy(IS1, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call ISDestroy(IS2, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  deallocate(dXsIndices)

  if (myid == 0) then 
     print *, 'Finished Mesh Initialization.'
  end if
  initializationSet = 1

end subroutine initializeWarping

subroutine getLoopFact(nLoop, fact)
  use constants
  use gridInput
  implicit none
  integer(kind=intType), intent(out) :: nLoop
  real(kind=realType), intent(out), dimension(3, 2) :: fact

  nLoop = 1
  fact = one
  if (iSymm > 0) then 
     nLoop = 2
     fact(iSymm, 2) = -one
  end if

end subroutine getLoopFact

