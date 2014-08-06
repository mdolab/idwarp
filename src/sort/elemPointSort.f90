
! Use the same number of points throughout
subroutine orderElemMerge(nPts,sortPoints,center,normal, a, middle, b)
  use precision
  implicit none

  ! subroutine variables
  integer(kind=intType):: nPts
  integer(kind=intType):: a,b,middle
  real(kind=realType),dimension(nPts,3):: sortPoints
  real(kind=realType),dimension(3)::center,normal
  
  ! Local variables
  integer(kind=intType):: ai,bi,ti,x,ierr
  real(kind=realType),dimension(nPts,3):: tmp
  real(kind=realType),dimension(2,3)::radialVec
  real(kind=realType),dimension(3)::cross
 
  ai = a
  bi = middle
  ti = a

  do while ((ai .lt. middle) .or. (bi .lt. b))
     radialVec(1,:) = sortPoints(ai+1,:)-center
     radialVec(2,:) = sortPoints(bi+1,:)-center
     call cross_product_3d(radialVec(1,:),radialVec(2,:),cross)
     if (ai .eq. middle) then
        tmp(ti+1,:) = sortPoints(bi+1,:)

        bi = bi + 1
     else if (bi .eq. b) then
        tmp(ti+1,:) = sortPoints(ai+1,:)
        ai = ai + 1
     ! else if (sortPoints(ai+1)%globalIndex .lt. sortPoints(bi+1)%globalIndex) then
     else if(dot_product(normal,cross).lt. 0)then
        tmp(ti+1,:) = sortPoints(ai+1,:)
        ai = ai + 1
     else
        tmp(ti+1,:) = sortPoints(bi+1,:)

        bi = bi + 1
     end if
     ti = ti + 1
  end do
  do x = a, b - 1
     sortPoints(x + 1,:) = tmp(x + 1,:)
  end do

end subroutine orderElemMerge

recursive subroutine orderElemPoints(nPts,points,center,normal, a, b)
  use precision
  implicit none
  integer(kind=intType):: nPts,a,b,diff
  real(kind=realType),dimension(nPts,3):: points
  real(kind=realType),dimension(3)::normal,center

  ! Local Variables
  
  diff = b - a
  !print *,'in orderElemPoints',diff,b,a
  if (diff .lt. 2) then
     return
  else
     diff = diff / 2
     !print *,'a',a,a+diff,diff,b
     call orderElemPoints(nPts,points,center,normal, a, a + diff)
     !nPtsNew = b-(a+diff)+1
     call orderElemPoints(nPts,points,center,normal, a + diff, b)
     call orderElemMerge(nPts,points,center,normal, a, a + diff, b)
     
  endif
end subroutine orderElemPoints

! subroutine indexShow(nPts,lst)
!   use precision
!   use gridData
!   implicit none

!   integer(kind=intType)::nPts
!   type(surfacePointType),dimension(nPts):: lst
!   integer(kind=intType):: x

!   do x = 1, nPts
!      print 100, 'index: ',x,' global: ',lst(x)%globalIndex,' elem: ',lst(x)%connectedElements(1,1),lst(x)%connectedElements(1,2),'loc',lst(x)%loc(1),lst(x)%loc(2),lst(x)%loc(3)
!   end do

! 100 format (A,i10,A,i10,A,2i10,A,3f10.6)
! end subroutine indexShow

