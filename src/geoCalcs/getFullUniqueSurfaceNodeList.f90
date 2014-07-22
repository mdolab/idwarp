! ====================================================================
! File: getFullUniqueSurfaceNodeList.f90
! Author: C.A.(Sandy) Mader
! Date Started: July 15, 2014
! Date Modified:


subroutine getFullUniqueSurfaceNodeList()
  use precision
  use gridData
  use communication
  implicit none
  include 'cgnslib_f.h'
  ! Local Variable
  integer(kind=intType)::zone,sec,ierr,idx
  integer(kind=intType):: surfSecCounter
  integer(kind=intType):: elem,conn,nConn,nElem, nPts
  type(surfacePointType),dimension(:),allocatable::tempSurfPoints, tempBoundaryPoints
  !type(surfacePointType),dimension(10)::testPoints
  integer(kind=intType)::nUniquePoints,surfCounter, boundaryCounter
  
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
        print *,'section',sec
        if(.not. gridDoms(zone)%isVolumeSection(sec))then
           nconn = gridDoms(zone)%surfaceSections(surfSecCounter)%nConn
           nElem = gridDoms(zone)%surfaceSections(surfSecCounter)%nElem
           if(gridDoms(zone)%surfaceSections(surfSecCounter)%isWallBC ) then
              nSurfNodes = nSurfNodes + nConn*nElem
           elseif(gridDoms(zone)%surfaceSections(surfSecCounter)%isBoundaryBC)then
              nBoundaryNodes = nBoundaryNodes + nConn*nElem
           end if
           surfSecCounter = surfSecCounter +1
        end if
     end do
  end do
           
  allocate(tempSurfPoints(nSurfNodes),tempBoundaryPoints(nBoundaryNodes),STAT=ierr)
  surfCounter=1
  boundaryCounter=1
  do zone = 1,nZones
     surfSecCounter = 1 
     do sec = 1,gridDoms(zone)%nSections
        print *,'section',sec
        if(.not. gridDoms(zone)%isVolumeSection(sec))then
           nconn = gridDoms(zone)%surfaceSections(surfSecCounter)%nConn
           nElem = gridDoms(zone)%surfaceSections(surfSecCounter)%nElem
           ! Loop over the surface and store the node data in the temp array
           do elem = 1,nElem
              do conn = 1,nConn
                 if(gridDoms(zone)%surfaceSections(surfSecCounter)%isWallBC ) then
                    tempSurfPoints(surfCounter)%globalIndex = &
                         gridDoms(zone)%surfaceSections(surfSecCounter)%elements(conn,elem)
                    tempSurfPoints(surfCounter)%proc = myID
                    tempSurfPoints(surfCounter)%zone = zone
                    tempSurfPoints(surfCounter)%surfSection = surfSecCounter
                    tempSurfPoints(surfCounter)%connectedElements(:) = -1_intType
                    tempSurfPoints(surfCounter)%connectedElements(1)=elem 
                    tempSurfPoints(surfCounter)%loc = &
                         gridDoms(zone)%points(tempSurfPoints(surfCounter)%globalIndex,:)
                    tempSurfPoints(surfCounter)%loc0 = &
                         gridDoms(zone)%points(tempSurfPoints(surfCounter)%globalIndex,:)
                    surfCounter = surfCounter + 1
                 elseif(gridDoms(zone)%surfaceSections(surfSecCounter)%isBoundaryBC)then
                    tempBoundaryPoints(boundaryCounter)%globalIndex = &
                         gridDoms(zone)%surfaceSections(surfSecCounter)%elements(conn,elem)
                    tempBoundaryPoints(boundaryCounter)%proc = myID
                    tempBoundaryPoints(boundaryCounter)%zone = zone
                    tempBoundaryPoints(boundaryCounter)%surfSection = surfSecCounter
                    tempBoundaryPoints(boundaryCounter)%connectedElements(:) = -1_intType
                    tempBoundaryPoints(boundaryCounter)%connectedElements(1)=elem 
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

  ! Sort the point array
  ! print *,'before sort',nConn,nElem,nPts
  ! call indexShow(nPts,tempPoints)
  ! print *,'sorting',nPts,shape(tempPoints)
  call pointMergeSort(nSurfNodes,tempSurfPoints,0_intType,nSurfNodes)
  call pointMergeSort(nBoundaryNodes,tempBoundaryPoints,0_intType,nBoundaryNodes)
  ! print *,'after sort'
  ! call indexShow(nPts,tempPoints)

  !compact the list
  call compactSurfacePointList(nSurfNodes,tempSurfPoints,nUniqueSurfPoints)
  call compactSurfacePointList(nBoundaryNodes,tempBoundaryPoints,nUniqueBoundaryPoints)
  print *,'Compacted sizes',nUniqueSurfPoints,nUniqueBoundaryPoints
  ! Allocate the permanent storage
  allocate(uniqueSurfaceNodes(nUniqueSurfPoints),STAT=ierr)
  allocate(uniqueBoundaryNodes(nUniqueBoundaryPoints),STAT=ierr)           
  uniqueSurfaceNodes = tempSurfPoints(1:nUniqueSurfPoints)
  uniqueBoundaryNodes = tempSurfPoints(1:nUniqueBoundaryPoints)
  
  deallocate(tempSurfPoints,tempBoundaryPoints,STAT=ierr)
  ! print *,'final points'
  ! call indexShow(nUniquePoints,gridDoms(zone)%surfaceSections(surfSecCounter)%uniqueSurfaceNodes)

end subroutine getFullUniqueSurfaceNodeList
