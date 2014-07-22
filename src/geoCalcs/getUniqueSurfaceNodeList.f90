! ====================================================================
! File: getUniqueSurfaceNodeList.f90
! Author: C.A.(Sandy) Mader
! Date Started: July 15, 2014
! Date Modified:


subroutine getUniqueSurfaceNodeList()
  use precision
  use gridData
  implicit none

  ! Local Variable
  integer(kind=intType)::zone,sec,ierr,idx
  integer(kind=intType):: surfSecCounter
  integer(kind=intType):: elem,conn,nConn,nElem, nPts
  type(surfacePointType),dimension(:),allocatable::tempPoints
  type(surfacePointType),dimension(10)::testPoints
  integer(kind=intType)::nUniquePoints
  ! begin execution
  
  ! loop over the surface sections
  ! Allocate a temporary work array the size of nConn*nElem
  ! Loop over the connectivity storing only unique nodes counting along the 
  ! way.
  ! once we have an accurate count, allocate a permanent array and store
  
  do zone = 1,nZones
     surfSecCounter = 1 
     do sec = 1,gridDoms(zone)%nSections
        print *,'section',sec
        if(.not. gridDoms(zone)%isVolumeSection(sec))then
           nconn = gridDoms(zone)%surfaceSections(surfSecCounter)%nConn
           nElem = gridDoms(zone)%surfaceSections(surfSecCounter)%nElem
           nPts = nConn*nElem
           !allocate(gridDoms(zone)%surfaceSections(surfSecCounter)%tmp(nConn*nElem),STAT=ierr
           allocate(tempPoints(nConn*nElem),STAT=ierr)
           ! Loop over the surface and store the node data in the temp array
           do elem = 1,nElem
              do conn = 1,nConn
                 idx = ((elem-1)*nConn)+conn
                 tempPoints(idx)%globalIndex = gridDoms(zone)%surfaceSections(surfSecCounter)%elements(conn,elem)
                 tempPoints(idx)%connectedElements(:) = -1_intType
                 tempPoints(idx)%connectedElements(1)=elem 
                 !tempPoints(idx)%loc = gridDoms(zone)%points(tempPoints(idx)%globalIndex,:)
              end do
           end do

           ! Sort the point array
           ! print *,'before sort',nConn,nElem,nPts
           ! call indexShow(nPts,tempPoints)
           ! print *,'sorting',nPts,shape(tempPoints)
           call pointMergeSort(nPts,tempPoints,0_intType,nPts)
           ! print *,'after sort'
           ! call indexShow(nPts,tempPoints)

           !compact the list
           call compactSurfacePointList(nPts,tempPoints,nUniquePoints)

           ! Allocate the permanent storage
           allocate(gridDoms(zone)%surfaceSections(surfSecCounter)%uniqueSurfaceNodes(nUniquePoints),STAT=ierr)
           
           gridDoms(zone)%surfaceSections(surfSecCounter)%nSurf = nUniquePoints
           
           gridDoms(zone)%surfaceSections(surfSecCounter)%uniqueSurfaceNodes = &
                tempPoints(1:nUniquePoints)

           deallocate(tempPoints,STAT=ierr)
           ! print *,'final points'
           ! call indexShow(nUniquePoints,gridDoms(zone)%surfaceSections(surfSecCounter)%uniqueSurfaceNodes)

           surfSecCounter = surfSecCounter +1
        end if
     end do
  end do
end subroutine getUniqueSurfaceNodeList
