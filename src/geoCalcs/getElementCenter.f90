! ====================================================================
! File: getElementCenter.f90
! Author: C.A.(Sandy) Mader
! Date Started: July 14, 2014
! Date Modified:

subroutine getElementCenter(nPts,nDim,points,center)
  use precision
  implicit none

  ! subroutine Variables
  integer(kind=intType)::nPts,nDim
  real(kind=realType),dimension(nPts,nDim)::points
  real(kind=realType),dimension(nDim)::center

  ! Local variables
  integer(kind=intType):: i,j
  real(kind=realType),dimension(nDim):: pointSum

  ! begin execution
  pointSum = 0
  do j = 1,nDim 
     do i = 1,nPts
        pointSum(j) = pointSum(j)+points(i,j)
     end do
     center(j) = pointSum(j)/nPts
  end do
  
end subroutine getElementCenter
