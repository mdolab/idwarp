! ====================================================================
! File: getFullUniqueSurfaceNodeList.f90
! Author: C.A.(Sandy) Mader
! Date Started: July 15, 2014
! Date Modified:


subroutine getFullUniqueSurfaceNodeList()
  use precision
  use gridData
  use sortData
  use communication
  implicit none
  include 'cgnslib_f.h'
  ! Local Variable
  integer(kind=intType)::zone,sec,ierr,idx
  integer(kind=intType):: surfSecCounter
  integer(kind=intType):: elem,conn,nConn,nElem, nPts,pt,globalIndex
  type(surfacePointType),dimension(:),allocatable::tempSurfPoints, tempBoundaryPoints
  type(surfacePointType),dimension(10)::testPoints
  integer(kind=intType)::nUniquePoints,surfCounter, boundaryCounter,nNodesTest
  integer(kind=intType)::elemIdx
  ! begin execution
  
  ! loop over all of the zones and sections.
  ! count up the total number of surface nodes and allocate an array of that size
  ! loop over the zones and sections again storing the informatoin on the nodes
  ! Run the compaction scheme to eliminate duplicates, then allocate the permanenet
  ! array and store

  nSurfNodes = 0
  nBoundaryNodes = 0
  do zone = 1,nZones
     surfSecCounter = 1 
     do sec = 1,gridDoms(zone)%nSections
        !print *,'surfsection',sec,gridDoms(zone)%isVolumeSection(sec)
        if(.not. gridDoms(zone)%isVolumeSection(sec))then
           !nconn = gridDoms(zone)%surfaceSections(surfSecCounter)%nConn
           nElem = gridDoms(zone)%surfaceSections(surfSecCounter)%nElem
           !print *,'nElem',nElem
           do elem = 1,nElem
              nConn = gridDoms(zone)%surfaceSections(surfSecCounter)%elemPtr(elem+1)-gridDoms(zone)%surfaceSections(surfSecCounter)%elemPtr(elem)

              if(gridDoms(zone)%surfaceSections(surfSecCounter)%isWallBC ) then
                 nSurfNodes = nSurfNodes + nConn
              elseif(gridDoms(zone)%surfaceSections(surfSecCounter)%isBoundaryBC)then
                 nBoundaryNodes = nBoundaryNodes + nConn
              end if
           end do
           surfSecCounter = surfSecCounter +1
        end if
     end do
  end do
  !print *,'nSurfNodes',nSurfNodes,nBoundaryNodes
!  stop
  allocate(tempSurfPoints(nSurfNodes),tempBoundaryPoints(nBoundaryNodes),STAT=ierr)
  surfCounter=1
  boundaryCounter=1
  do zone = 1,nZones
     surfSecCounter = 1 
     do sec = 1,gridDoms(zone)%nSections
        !print *,'tempPointSection',sec,surfSecCounter,gridDoms(zone)%isVolumeSection(sec)
        if(.not. gridDoms(zone)%isVolumeSection(sec))then
           !nconn = gridDoms(zone)%surfaceSections(surfSecCounter)%nConn
           nElem = gridDoms(zone)%surfaceSections(surfSecCounter)%nElem
           ! Loop over the surface and store the node data in the temp array
           !print *,'nElem',nElem
           !print *,'Boundary Type',gridDoms(zone)%surfaceSections(surfSecCounter)%isWallBC,gridDoms(zone)%surfaceSections(surfSecCounter)%isBoundaryBC,gridDoms(zone)%surfaceSections(surfSecCounter)%isSymmBC
           do elem = 1,nElem
              nConn = gridDoms(zone)%surfaceSections(surfSecCounter)%elemPtr(elem+1)-gridDoms(zone)%surfaceSections(surfSecCounter)%elemPtr(elem)
              !print *,'nConn',nConn
              do conn = 1,nConn
                 if(gridDoms(zone)%surfaceSections(surfSecCounter)%isWallBC ) then
                    elemIdx = gridDoms(zone)%surfaceSections(surfSecCounter)%elemPtr(elem)
                    tempSurfPoints(surfCounter)%globalIndex = &
                         gridDoms(zone)%surfaceSections(surfSecCounter)%elemConn(elemIdx+conn-1)
                    tempSurfPoints(surfCounter)%proc = myID
                    tempSurfPoints(surfCounter)%zone = zone
                    !tempSurfPoints(surfCounter)%surfSection = surfSecCounter
                    tempSurfPoints(surfCounter)%connectedElements(:,:) = -1_intType
                    tempSurfPoints(surfCounter)%connectedElements(1,1)=elem 
                    tempSurfPoints(surfCounter)%connectedElements(1,2)=surfSecCounter
                    !print *,'surfCounter',surfCounter,tempSurfPoints(surfCounter)%connectedElements(1,2),tempSurfPoints(surfCounter)%globalIndex,elemIdx+conn-1
                    tempSurfPoints(surfCounter)%loc = &
                         gridDoms(zone)%points(tempSurfPoints(surfCounter)%globalIndex,:)
                    tempSurfPoints(surfCounter)%loc0 = &
                         gridDoms(zone)%points(tempSurfPoints(surfCounter)%globalIndex,:)
                    surfCounter = surfCounter + 1
                 elseif(gridDoms(zone)%surfaceSections(surfSecCounter)%isBoundaryBC)then
                    !print *,'elem',elem,shape(gridDoms(zone)%surfaceSections(surfSecCounter)%elemPtr)
                    elemIdx = gridDoms(zone)%surfaceSections(surfSecCounter)%elemPtr(elem)
                    !print *,'elemidx',elemIdx,conn,nElem,shape(gridDoms(zone)%surfaceSections(surfSecCounter)%elemConn)
                    !print *,'boundary Counter',boundaryCounter,shape(tempBoundaryPoints),nBoundaryNodes
                    tempBoundaryPoints(boundaryCounter)%globalIndex = &
                         gridDoms(zone)%surfaceSections(surfSecCounter)%elemConn(elemIdx+conn-1)
                    tempBoundaryPoints(boundaryCounter)%proc = myID
                    tempBoundaryPoints(boundaryCounter)%zone = zone
                    !tempBoundaryPoints(boundaryCounter)%surfSection = surfSecCounter
                    tempBoundaryPoints(boundaryCounter)%connectedElements(:,:) = -1_intType
                    tempBoundaryPoints(boundaryCounter)%connectedElements(1,1)=elem 
                    tempBoundaryPoints(boundaryCounter)%connectedElements(1,2)=surfSecCounter
                    !print *,'boundaryPoint',boundaryCounter,tempBoundaryPoints(boundaryCounter)%connectedElements(1,2)
                    tempBoundaryPoints(boundaryCounter)%loc = &
                         gridDoms(zone)%points(tempBoundaryPoints(boundaryCounter)%globalIndex,:)
                    tempBoundaryPoints(boundaryCounter)%loc0 = &
                         gridDoms(zone)%points(tempBoundaryPoints(boundaryCounter)%globalIndex,:)
                    boundaryCounter = boundaryCounter + 1
                 end if
                 ! !idx = ((elem-1)*nConn)+conn
                 ! tempPoints(idx)%globalIndex = gridDoms(zone)%surfaceSections(surfSecCounter)%elements(conn,elem)
                 ! tempPoints(idx)%connectedElements(:) = -1_intType
                 ! tempPoints(idx)%connectedElements(1)=elem 
                 ! !tempPoints(idx)%loc = gridDoms(zone)%points(tempPoints(idx)%globalIndex,:)
              end do
           end do
           surfSecCounter = surfSecCounter +1
        end if
     end do
  end do
  !print *,'boundary Counter',boundaryCounter,surfCounter
  

  ! ! Check that the global indices match up
  ! surfCounter=1
  ! do zone = 1,nZones
  !    do pt = 1, 100!nSurfNodes
  !       globalIndex = tempSurfPoints(pt)%globalIndex
  !       print *,pt,'surfpoints',tempSurfPoints(pt)%loc,gridDoms(zone)%points(globalIndex,:),globalIndex
  !    end do
  ! end do
  ! !stop
  ! ! Sort the point array
  ! nNodesTest = 10
  ! testPoints = tempSurfPoints(1:10)
  ! call indexShow(nNodesTest,testPoints)
  ! ! print *,'sorting',nPts,shape(tempPoints)
  ! allocate(sortPoints(nNodesTest),STAT=ierr)
  ! sortPoints = testPoints
  ! call pointMergeSort(nNodesTest,0_intType,nNodesTest)
  ! testPoints = sortPoints
  ! deallocate(sortPoints,STAT=ierr)
  ! print *,'after sort'
  ! call indexShow(nNodesTest,testPoints)
  ! stop

  ! print *,'before sort',nConn,nElem,nSurfNodes
  ! call indexShow(nSurfNodes,tempSurfPoints)
  ! print *,'sorting',nPts,shape(tempPoints)
  allocate(sortPoints(nSurfNodes),STAT=ierr)
  sortPoints = tempSurfPoints
  !call pointMergeSort(nSurfNodes,tempSurfPoints,0_intType,nSurfNodes)
  call pointMergeSort(nSurfNodes,0_intType,nSurfNodes)
  tempSurfPoints = sortPoints
  deallocate(sortPoints,STAT=ierr)
  ! print *,'after sort'
  ! call indexShow(nSurfNodes,tempSurfPoints)

  ! repeat for Boundary
  allocate(sortPoints(nBoundaryNodes),STAT=ierr)
  sortPoints = tempBoundaryPoints
  !call pointMergeSort(nBoundaryNodes,tempBoundaryPoints,0_intType,nBoundaryNodes)
  call pointMergeSort(nBoundaryNodes,0_intType,nBoundaryNodes)
  tempBoundaryPoints = sortPoints
  deallocate(sortPoints,STAT=ierr)

  ! do zone = 1,nZones
  !    do pt = 1, 100!nSurfNodes
  !       globalIndex = tempSurfPoints(pt)%globalIndex
  !       print *,pt,'Sortedsurfpoints',tempSurfPoints(pt)%loc,gridDoms(zone)%points(globalIndex,:),globalIndex
  !    end do
  ! end do

  !compact the list
  call compactSurfacePointList(nSurfNodes,tempSurfPoints,nUniqueSurfPoints)
  call compactSurfacePointList(nBoundaryNodes,tempBoundaryPoints,nUniqueBoundaryPoints)
  !print *,'Compacted sizes',nUniqueSurfPoints,nUniqueBoundaryPoints
  ! Allocate the permanent storage
  allocate(uniqueSurfaceNodes(nUniqueSurfPoints),STAT=ierr)
  allocate(uniqueBoundaryNodes(nUniqueBoundaryPoints),STAT=ierr)           
  uniqueSurfaceNodes = tempSurfPoints(1:nUniqueSurfPoints)
  uniqueBoundaryNodes = tempBoundaryPoints(1:nUniqueBoundaryPoints)
  
  deallocate(tempSurfPoints,tempBoundaryPoints,STAT=ierr)
  ! print *,'final points'
  ! call indexShow(nUniqueSurfPoints,uniqueSurfaceNodes)
  ! stop
  ! Now loop over the nodes and mark the boundary nodes in the global point list
  ! First loop over the Farfield boundaries
  do pt = 1, nUniqueBoundaryPoints
     zone = uniqueBoundaryNodes(pt)%zone
     globalIndex = uniqueBoundaryNodes(pt)%globalIndex
     gridDoms(zone)%isSurfaceNode(globalIndex) = .True.      
  end do
     ! Then loop over the surface boundaries
  do pt = 1, nUniqueSurfPoints
     zone = uniqueSurfaceNodes(pt)%zone
     globalIndex = uniqueSurfaceNodes(pt)%globalIndex
     gridDoms(zone)%isSurfaceNode(globalIndex) = .True.  
  end do
  ! do pt = 1, 100!nUniqueSurfPoints
  !    zone = uniqueSurfaceNodes(pt)%zone
  !    globalIndex = uniqueSurfaceNodes(pt)%globalIndex
  !    print *,'surf',pt,globalIndex,gridDoms(zone)%points(globalIndex,:),uniqueSurfaceNodes(pt)%loc
  !    !gridDoms(zone)%points(globalIndex,:) = uniqueSurfaceNodes(pt)%loc  
  ! end do
  ! stop
end subroutine getFullUniqueSurfaceNodeList
