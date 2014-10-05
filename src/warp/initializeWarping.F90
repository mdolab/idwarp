subroutine initializeWarping(pts, ndof, faceSizesLocal, faceConnLocal, nFaceSizesLocal, &
     nFaceConnLocal)
  
  use gridData
  use gridInput
  use communication
  use kd_tree

  implicit none
 
  ! Input
  real(kind=realType), intent(in) :: pts(ndof)
  integer(kind=intType), intent(in) :: ndof
  integer(kind=intType), dimension(nFaceSizesLocal), intent(in) :: faceSizesLocal
  integer(kind=intTYpe), dimension(nFaceConnLocal), intent(in) :: faceConnLocal
  integer(kind=intType), intent(in) :: nFaceSizesLocal, nFaceConnLocal

  ! Working
  integer(kind=intType) :: nNodesTotal, ierr, iStart, iEnd, iProc
  integer(kind=intType) :: i, j, k, n, ii, nPts, ind, curInd, maxConnectedFace
  integer(kind=intType), dimension(:), allocatable :: nNodesProc, cumNodesProc
  real(kind=realType),dimension(:,:), allocatable :: allNodes
  integer(kind=intType), dimension(:), allocatable :: faceSizesMirror, faceConnMirror
  real(kind=realType), dimension(:, :), allocatable :: allMirrorPts, uniquePts
  real(kind=realType), dimension(:), allocatable :: costs, procEndCosts
  integer(kind=intType), dimension(:), allocatable :: surfSizesProc, surfSizesDisp
  integer(kind=intType), dimension(:), allocatable :: surfConnProc, surfConnDisp
  integer(kind=intType), dimension(:), allocatable :: link, tempInt, invIndices
  integer(kind=intType), dimension(:), allocatable :: procSplits, procSplitsLocal
  real(kind=realType), dimension(:), allocatable :: denomenator0Copy
  integer(kind=intType), Pointer :: indices(:)
  real(kind=realtype) :: fact(3), xcen(3), dx(3), r(3), costOffset
  real(Kind=realType) :: totalCost, averageCost
  integer(kind=intType) :: nFaceConnMirror, nFaceSizesMirror, nMirrorNodes, nVol
  integer(kind=intType) :: curBin, newBin, newDOFProc, totalDOF

   interface 
     subroutine pointReduce(pts, N, tol, uniquePts, link, nUnique)
       use precision
       implicit none
       real(kind=realType), intent(in), dimension(:, :) :: pts
       integer(kind=intType), intent(in) :: N
       real(kind=realType), intent(in) :: tol
       real(kind=realType), intent(out), dimension(:, :) :: uniquePts
       integer(kind=intType), intent(out), dimension(:) :: link
       integer(kind=intType), intent(out) :: nUnique
     end subroutine pointReduce
  end interface

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

  ! ----------------------------------------------------------------------
  !   Step 3: Deal with the symmetry plane
  ! ----------------------------------------------------------------------
  if (iSymm > 0) then 

     ! Mirror the nodes according to iSymm
     allocate(allMirrorPts(3, 2*nNodesTotal))
     fact = one
     fact(iSymm) = -one

     ! Copy the nodes
     j = nNodesTotal
     do i=1, nNodesTotal
        allMirrorPts(:, i) = allNodes(:, i)
        allMirrorPts(:, i + j) = allNODes(:, i)*fact
     end do

     ! And now we do the sizes/connectivity
     allocate(faceSizesMirror(2*nFaceSizesLocal), &
              faceConnMirror(2*nFaceConnLocal))
     ! Face sizes is a direct copy
     do i=1, nFaceSizesLocal
        faceSizesMirror(i) = faceSizesLocal(i)
        faceSizesMirror(i+nFaceSizesLocal) = faceSizesLocal(i)
     end do

     ! Connectivity is a little tricker since we need need to switch
     ! the orientation of the face nodes to keep everything consistent
     k = 0
     n = nFaceConnLocal
     do i=1, nFaceSizesLocal
        nFace = faceSizesLocal(i)
        do j=1, nFace
           faceConnMirror(k + j) = faceConnLocal(k + j)
           faceConnMirror(k + n + nFace - j + 1) = faceConnLocal(k + j) + nNodesTotal
        end do
        k = k + nFace
     end do
     nFaceSizesMirror = nFaceSizesLocal*2
     nFaceConnMirror = nFaceConnLocal*2
     nMirrorNodes = nNodesTotal*2
  else
     ! Just copy in the values from all Nodes, the connectivity is
     ! already done
     allocate(allMirrorPts(3, nNodesTotal))
     do i=1, nNodesTotal
        allMirrorPts(:, i) = allNodes(:, i)
     end do

     allocate(faceSizesMirror(nFaceSizesLocal), &
              faceConnMirror(nFaceConnLocal))
     faceSizesMirror = faceSizesLocal
     faceConnMirror = faceConnLocal
     nFaceSizesMirror = nFaceSizesLocal
     nFaceConnMirror = nFaceConnLocal
     nMirrorNodes = nNodesTotal
  end if
  deallocate(allNodes)

  ! ----------------------------------------------------------------------
  !   Step 4: Gather up the connectivity information
  ! ----------------------------------------------------------------------

  ! Gather the displacements
  allocate(surfSizesProc(nProc), surfSizesDisp(0:nProc), &
           surfConnProc(nProc), surfConnDisp(0:nProc))

  call MPI_allgather(nFaceSizesMirror, 1, MPI_INTEGER, surfSizesProc, 1, MPI_INTEGER, &
       warp_comm_world, ierr)
  call EChk(ierr,__FILE__,__LINE__)
  
  call MPI_allgather(nFaceConnMirror, 1, MPI_INTEGER, surfConnProc, 1, MPI_INTEGER, &
       warp_comm_world, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Finish the displacment calc:
  surfSizesDisp(:) = 0
  surfConnDisp(:) = 0
  do i=1, nProc
     surfSizesDisp(i) = surfSizesDisp(i-1) + surfSizesProc(i)
     surfConnDisp(i) = surfConnDisp(i-1) + surfConnProc(i)
  end do
  nFace = surfSizesDisp(nProc)
  lenFaceConn = surfConnDisp(nProc)
  allocate(facePtr(0:nFace), faceConn(lenFaceConn))

  ! And now for the super magical all-gather-v's!
    call MPI_Allgatherv(&
       faceSizesMirror, nFaceSizesMirror, MPI_INTEGER, &
       facePtr(1:)   , surfSizesProc, surfSizesDisp, MPI_INTEGER, warp_comm_world, ierr)
  call EChk(ierr,__FILE__,__LINE__)

    call MPI_Allgatherv(&
         faceConnMirror, nFaceConnMirror, MPI_INTEGER, &
         faceConn     , surfConnProc, surfConnDisp, MPI_INTEGER, warp_comm_world, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Note that we put the result of the "faceSizes" allgather into
  ! 'facePtr' which implies it should be size pointer...which we will
  ! make it into now. Note that facePtr is 0-based which makes things
  ! easier. 
  facePtr(0) = 0
  do i=2, nFace
     facePtr(i) = facePtr(i) + facePtr(i-1)
  end do

  ! We now have to update faceConn since it is currently uses the
  ! indexing from the mirror coordinates so we need up update by the
  ! offset in the Xs nodes
  do iProc=0, nProc-1
     do i=surfConnDisp(iProc)+1, surfConnDisp(iProc+1)
        faceConn(i) = faceConn(i) + cumNodesProc(iProc)
     end do
  end do

  allocate(uniquePts(3, nMirrorNodes), link(nMirrorNodes))
  call pointReduce(allMirrorPts, nMirrorNodes, symmTol, uniquePts, link, nUnique)

  ! Let the user know the actual number of Surface nodes used for the calc
  if (myid == 0) then
     write(*,"(a)", advance="no") '#--------------------------------#'
     print "(1x)"  
     write(*,"(a)", advance="no") " Unique Surface Nodes : "
     write(*,"(I9,1x)",advance="no") nUnique
     print "(1x)"  
     write(*,"(a)", advance="no") '#--------------------------------#'
     print "(1x)"   
  end if

  ! Before we create the final Xu and friends we will do one final
  ! aditional mapping: We will create the KD to get the indices that
  ! *would* sort the tree. By doing this first, this removes the
  ! additional layer of indirection during the tree traversals of the
  ! nodal properties to be interpolated. Since this will be done
  ! BILLIONS and BILLIONS of times, it is worth making these accesses
  ! fast and essentially in order in memory.
  mytree=>create_tree(uniquePts(:, 1:nUnique))
  call labelNodes(mytree)
  call setCenters(mytree)

  indices => mytree%indexes
  ! Compute the inverse of the indices
  allocate(invIndices(nUnique))
  do i=1, nUnique
     invIndices(indices(i)) = i
  end do

  ! Now we can go through and update the conn. Note that conn is
  ! currently 0-based so we need the plus one when indexing into link
  do i=1, lenFaceConn
     faceConn(i) = invIndices(link(faceConn(i) + 1))
  end do

  ! We also need to compute the nodeToElem pointer....This data
  ! structure contains the number of elements each node is connected
  ! to as well as the indices of the faces. It is stored as 2D array
  ! for simplicitity...we first have to determine the maximum number
  ! of cells connected to each node.

  allocate(tempInt(nUnique))
  tempInt = 0
  do i=1, nFace
     nPts = facePtr(i) - facePtr(i-1)
  
     do j=1, nPts
        ind = faceConn(facePtr(i-1) + j)
        tempInt(ind) = tempInt(ind) + 1
     end do
  end do

  maxConnectedFace = maxval(tempInt)
  allocate(nodeToElem(maxConnectedFace + 1, nUnique))
  deallocate(tempInt)

  ! Now we loop back again and actually assign the faces to the elements
  nodeToElem = 0
  do i=1, nFace
     nPts = facePtr(i) - facePtr(i-1)
     do j=1, nPts
        ind = faceConn(facePtr(i-1) + j)
        nodeToElem(1, ind) = nodeToElem(1, ind) + 1
        curInd = nodeToElem(1, ind)
        nodeToElem(curInd+1, ind) = i
     end do
  end do

  ! Now we need to compute the one-time mapping that goes between the
  ! unique set of surface nodes and Xs. This will be required for the
  ! sensitivity calc. 
  allocate(Xu(3, nUnique), Xu0(3, nUnique), XuInd(nUnique), XuFact(3, nUnique))

  ! Copy uniquePoints into Xu0 which will be fixed for all time:
  do i=1, nUnique
     Xu0(:, i) = uniquePts(:, indices(i))
  end do

  ! Use zero to indicate an un-assigned node
  XuInd(:) = 0

  fact(:) = one
  if (iSymm > 1) then
     fact(iSymm) = -one
  end if

  do i=1, nMirrorNodes
     if (i <= nNodesTotal) then ! Regular part:
        if (XuInd(invindices(link(i))) == 0) then ! Un-assigned
           XuInd(invindices(link(i))) = i
           XuFact(:, invindices(link(i))) = one
        end if
     else ! Mirror part:
        if (XuInd(invindices(link(i))) == 0) then ! Un-assigned
           XuInd(invindices(link(i))) = i - nNodesTotal
           XuFact(:, invindices(link(i))) = fact
        end if
     end if
  end do

  ! Compute Ldef based on the size of the mesh
  Xcen = zero
  do i=1,nUnique
     Xcen = Xcen + Xu0(:, i)
  end do
  Xcen = Xcen / nUnique

  ! Now get the max dist from Xcen to any node:
  Ldef0 = zero
  do i=1,nUnique
     dx = Xu0(:, i) - Xcen
     Ldef0 = max(Ldef0, sqrt(dx(1)**2 + dx(2)**2 + dx(3)))
  end do

  ! Allocate the space for Mi, Bi, and Ai
  allocate(Mi(3, 3, nUnique), Bi(3, nUnique), Ai(nUnique))
  allocate(normals(3, nUnique), normals0(3, nUnique))
  
  ! Run the compute nodal properties to initialize the normal vectors.
  call initNodalProperties()

  ! Set the input parameters to the module
  mytree%Ldef = Ldef0 * LdefFact
  mytree%alphaToBexp = alpha**bExp
  mytree%errTol = errTol

  ! Perform some initialization on the tree
   call setAi(mytree, Ai)
   call setXu0(mytree, Xu0)
   call labelNodes(mytree)
   call computeErrors(mytree)

  ! Compute the approximate denomenator:
  call VecGetArrayF90(commonGridVec, Xv0Ptr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  if (myid == 0) then 
     print *, 'Computing Denomenator Estimate...'
  end if

  nVol = size(Xv0Ptr)/3
  allocate(denomenator0(nVol))

  ! Compute the (approx) denomenator. This needs to be done only once.
  wiLoop: do j=1, nVol
     call getWiEstimate(mytree, Xv0Ptr(3*j-2:3*j), denomenator0(j))
  end do wiLoop

  ! If we want to load balance, which generally should be done, we
  ! have to do a dry run loop that determines just the number of
  ! ops that a mesh warp *would* use. This integer is scaled by the
  ! "brute force" cost, so each individual cost will be less than
  ! 1.0
  if (myid == 0) then 
     print *, 'Load Balancing...'
  end if
  allocate(costs(nVol))
  dryRunLoop: do j=1, nVol
     call dryRun(mytree, Xv0Ptr(3*j-2:3*j), denomenator0(j), ii)
     costs(j) = dble(ii)/nUnique
  end do dryRunLoop
  
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
  ! gridVec. So we dump it into numerator, which *is* the right size
  ! and the copy after

 ! We can now also allocate the final space for the numerator and denomenator
  allocate(numerator(3, newDOFProc/3), denomenator(newDOFProc/3))

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
  
  ! Copy denomenator0 out of the numerator
  do i=1,newDOFProc/3
     denomenator0(i) = numerator(1, i)
  end do
  warpMeshDOF = newDOFProc

  ! Deallocate the memory from this subroutine
  deallocate(nNodesProc, cumNodesProc)
  deallocate(faceSizesMirror, faceConnMirror, allMirrorPts, uniquePts)
  deallocate(surfSizesProc, surfSizesDisp, surfConnProc, surfConnDisp)
  deallocate(link)
  nullify(indices)
  deallocate(invIndices)
  deallocate(procEndCosts, procSplits, procSplitsLocal, denomenator0Copy)
  if (myid == 0) then 
     print *, 'Finished Mesh Initialization.'
  end if

end subroutine initializeWarping

