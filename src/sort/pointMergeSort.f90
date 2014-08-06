! Use the same number of points throughout
!subroutine pointMerge(nPts,sortPoints, a, middle, b)
subroutine pointMerge(nPts, a, middle, b)
  use precision
  use gridData
  use sortData
  implicit none

  ! subroutine variables
  integer(kind=intType):: nPts
  integer(kind=intType):: a,b,middle
  !type(surfacePointType),dimension(nPts):: sortPoints
  
  ! Local variables
  integer(kind=intType):: ai,bi,ti,x,ierr
  !type(surfacePointType),dimension(nPts):: tmp
  type(surfacePointType),dimension(:), allocatable:: tmp
  
  allocate(tmp(nPts),STAT=ierr)
  
  ai = a
  bi = middle
  ti = a

  do while ((ai .lt. middle) .or. (bi .lt. b))
     if (ai .eq. middle) then
        tmp(ti+1) = sortPoints(bi+1)
        !call copySurfacePoint(sortPoints(bi+1),tmp(ti+1))
        bi = bi + 1
     else if (bi .eq. b) then
        tmp(ti+1) = sortPoints(ai+1)
        !call copySurfacePoint(sortPoints(ai+1),tmp(ti+1))
        ai = ai + 1
     else if (sortPoints(ai+1)%globalIndex .lt. sortPoints(bi+1)%globalIndex) then
        tmp(ti+1) = sortPoints(ai+1)
        !print *,'ai',ai+1,ti+1,shape(sortPoints),shape(tmp)
        !call copySurfacePoint(sortPoints(ai+1),tmp(ti+1))
        ai = ai + 1
     else
        tmp(ti+1) = sortPoints(bi+1)
        !call copySurfacePoint(sortPoints(bi+1),tmp(ti+1))
        bi = bi + 1
     end if
     ti = ti + 1
  end do
  do x = a, b - 1
     sortPoints(x + 1) = tmp(x + 1)
     !call copySurfacePoint(tmp(x+1),sortPoints(x+1))
  end do
  deallocate(tmp,STAT=ierr)
end subroutine pointMerge

!recursive subroutine pointMergeSort(nPts,points, a, b)
recursive subroutine pointMergeSort(nPts, a, b)
  use precision
  use gridData
  implicit none
  integer(kind=intType):: nPts,a,b,diff
  !type(surfacePointType),dimension(nPts)::points

  ! Local Variables
!  integer(kind=intType)::nPtsNew

  
  diff = b - a
  !print *,'in pointMergeSort',diff,b,a
  if (diff .lt. 2) then
     return
  else
     diff = diff / 2
     !print *,'a',a,a+diff,diff,b
     !call pointMergeSort(nPts,points, a, a + diff)
     call pointMergeSort(nPts, a, a + diff)
     !nPtsNew = b-(a+diff)+1
     !call pointMergeSort(nPts,points, a + diff, b)
     call pointMergeSort(nPts, a + diff, b)
     !call pointMerge(nPts,points, a, a + diff, b)
     call pointMerge(nPts, a, a + diff, b)
  endif
end subroutine pointMergeSort

subroutine indexShow(nPts,lst)
  use precision
  use gridData
  implicit none

  integer(kind=intType)::nPts
  type(surfacePointType),dimension(nPts):: lst
  integer(kind=intType):: x

  do x = 1, nPts
     print 100, 'index: ',x,' global: ',lst(x)%globalIndex,' elem: ',lst(x)%connectedElements(1,1),lst(x)%connectedElements(1,2),'loc',lst(x)%loc(1),lst(x)%loc(2),lst(x)%loc(3)
  end do

100 format (A,i10,A,i10,A,2i10,A,3f10.6)
end subroutine indexShow

! subroutine copySurfacePoint(ptIn,ptOut)
!   use precision
!   use gridData
!   implicit none

!   ! Subroutine Variables
!   type(surfacePointType), intent(in)::ptIn
!   type(surfacePointType), intent(out)::ptOut

!   !begin Execution
!   ptOut%globalIndex = ptIn%globalIndex
!   ptOut%proc = ptIn%proc
!   ptOut%zone = ptIn%zone
!   ptOut%connectedElements = ptIn%connectedElements
!   ptOut%loc = ptIn%loc
!   ptOut%loc0 = ptIn%loc0
!   ptOut%normal = ptIn%normal
!   ptOut%normal0 = ptIn%normal0
!   ptOut%Mi = ptIn%Mi
!   ptOut%bi = ptIn%bi
!   ptOut%Ai = ptIn%Ai
!   ptOut%Wi = ptIn%Wi
!   ptOut%Si = ptIn%Si

! end subroutine copySurfacePoint
