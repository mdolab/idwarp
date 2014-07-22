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
  integer(kind=intType)::i,surfSecCounter
  integer(kind=intType)::currentSlot,currentNodeNumber,currElem
  integer(kind=intType)::nodeNumber, connectedElement

  !begin execution

  ! Loop over the points in the list 
  ! keep a counter of the current node number
  ! if this matches the current node, copy the connected element number then erase
  ! if this is a new node, copy to the next available spot and move on.

   !Loop over points
  currentSlot = 1
  currentNodeNumber = 1
  currElem = 1
  do i = 1,nPts
     nodeNumber = tempPoints(i)%globalIndex
     connectedElement =  tempPoints(i)%connectedElements(1)
!     print *,'current Numbers',nodeNumber,connectedElement
     if (nodeNumber .eq. currentNodeNumber)then
        !print *,'eq'
        ! Current node is part of current set.
        tempPoints(currentSlot)%connectedElements(currElem) = connectedElement
        currElem = currElem + 1
     else
        !print *,'incrementing'
        ! We have moved to the next global node number, increment the current slot
        ! and reset the currElem counter
        currentSlot = currentSlot + 1
        currElem = 1
        
        ! We have already sorted the entries in tempPoints by global index so
        ! now set the current node number to matches the global index 
        ! of the current point and then proceed.
        currentNodeNumber = nodeNumber
        !Store all of the data
        tempPoints(currentSlot)%connectedElements(:) = -1_intType
        tempPoints(currentSlot)%connectedElements(currElem) = connectedElement
        !tempPoints(currentSlot)%loc = tempPoints(i)%loc
        tempPoints(currentSlot)%globalIndex = currentNodeNumber
        currElem = currElem + 1
     end if
     !print 100,'Updated',i,currentNodeNumber,currElem-1,' elems ',tempPoints(currentSlot)%connectedElements,nodeNumber,currentSlot
  end do
  
  nUniquePts = currentSlot
  
100 format (A,i7,i7,i7,A,10i5,i10,i10)
end subroutine compactSurfacePointList
