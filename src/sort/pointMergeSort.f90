! Use the same number of points throughout
subroutine pointMerge(nPts,lst, a, middle, b)
  use precision
  use gridData
  implicit none

  ! subroutine variables
  integer(kind=intType):: nPts
  integer(kind=intType):: a,b,middle
  type(surfacePointType),dimension(nPts):: lst
  
  ! Local variables
  integer(kind=intType):: ai,bi,ti,x
  type(surfacePointType),dimension(nPts):: tmp

  ai = a
  bi = middle
  ti = a

  do while ((ai .lt. middle) .or. (bi .lt. b))
     if (ai .eq. middle) then
        tmp(ti+1) = lst(bi+1)
        bi = bi + 1
     else if (bi .eq. b) then
        tmp(ti+1) = lst(ai+1)
        ai = ai + 1
     else if (lst(ai+1)%globalIndex .lt. lst(bi+1)%globalIndex) then
        tmp(ti+1) = lst(ai+1)
        ai = ai + 1
     else
        tmp(ti+1) = lst(bi+1)
        bi = bi + 1
     end if
     ti = ti + 1
  end do
  do x = a, b - 1
     lst(x + 1) = tmp(x + 1)
  end do

end subroutine pointMerge

recursive subroutine pointMergeSort(nPts,points, a, b)
  use precision
  use gridData
  implicit none
  integer(kind=intType):: nPts,a,b,diff
  type(surfacePointType),dimension(nPts)::points

  ! Local Variables
  integer(kind=intType)::nPtsNew

  diff = b - a

  if (diff .lt. 2) then
     return
  else
     diff = diff / 2
     call pointMergeSort(nPts,points, a, a + diff)
     nPtsNew = b-(a+diff)+1
     call pointMergeSort(nPts,points, a + diff, b)
     call pointMerge(nPts,points, a, a + diff, b)
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
     print 100, 'index: ',x,' global: ',lst(x)%globalIndex,' elem: ',lst(x)%connectedElements(1)!,lst(x)%loc
  end do

100 format (A,i10,A,i10,A,i10)!,3f8.2)
end subroutine indexShow
