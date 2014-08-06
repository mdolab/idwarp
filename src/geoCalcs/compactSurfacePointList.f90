! ====================================================================
! File: compactSurfacePointList.f90
! Author: C.A.(Sandy) Mader
! Date Started: July 18, 2014
! Date Modified:

subroutine compactSurfacePointList(nPts,tempPoints,nUniquePts)

  use precision
  use gridData
  implicit none

  !Subroutine Variables
  integer(kind=intType)::nPts
  integer(kind=intType),intent(out)::nUniquePts
  
  type(surfacePointType),dimension(nPts)::tempPoints

  ! Local Variables
  integer(kind=intType)::i,surfSecCounter,j
  integer(kind=intType)::currentSlot,currentNodeNumber,currElem
  integer(kind=intType)::nodeNumber
  integer(kind=intType),dimension(2)::connectedElement

  !begin execution

  ! Loop over the points in the list 
  ! keep a counter of the current node number
  ! if this matches the current node, copy the connected element number then erase
  ! if this is a new node, copy to the next available spot and move on.

   !Loop over points
  currentSlot = 1
  !print *,'tempPoints',shape(tempPoints)
  currentNodeNumber = tempPoints(1)%globalIndex!1
  currElem = 1
  do i = 1,nPts
     nodeNumber = tempPoints(i)%globalIndex
     connectedElement =  tempPoints(i)%connectedElements(1,:)
     ! if (connectedElement(2) > 6) then
     !    print *,'current Numbers',nodeNumber,connectedElement
     ! end if
     if (nodeNumber .eq. currentNodeNumber)then
        !print *,'eq',currentSlot,currElem,nodeNumber,currentNodeNumber
        ! Current node is part of current set.
        tempPoints(currentSlot)%connectedElements(currElem,:) = connectedElement
        currElem = currElem + 1
     else
        !print *,'incrementing',currentSlot
        ! We have moved to the next global node number, increment the current slot
        ! and reset the currElem counter
        currentSlot = currentSlot + 1
        currElem = 1
        
        ! We have already sorted the entries in tempPoints by global index so
        ! now set the current node number to matches the global index 
        ! of the current point and then proceed.
        currentNodeNumber = nodeNumber
        !Store all of the data
        tempPoints(currentSlot)%connectedElements(:,:) = -1_intType
        tempPoints(currentSlot)%connectedElements(currElem,:) = connectedElement
        tempPoints(currentSlot)%loc = tempPoints(i)%loc
        tempPoints(currentSlot)%loc0 = tempPoints(i)%loc0
        tempPoints(currentSlot)%globalIndex = currentNodeNumber
        currElem = currElem + 1
     end if
     !print 100,'Updated',i,currentNodeNumber,currElem-1,' elems ',tempPoints(currentSlot)%connectedElements,nodeNumber,currentSlot
  end do
  
  nUniquePts = currentSlot
!   do i = 1,nUniquePts
!      do j = 1,30
!         if (tempPoints(i)%connectedElements(j,1).ge.0)then
!            if(tempPoints(i)%connectedElements(j,2)>6) then
!               print *,'Compacted list',i,j,tempPoints(i)%connectedElements(j,:)
!            end if
!         end if
!      end do
!   end do
! stop
100 format (A,i7,i7,i7,A,10i5,i10,i10)
end subroutine compactSurfacePointList
