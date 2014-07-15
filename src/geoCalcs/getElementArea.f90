! ====================================================================
! File: getElementArea.f90
! Author: C.A.(Sandy) Mader
! Date Started: July 14, 2014
! Date Modified:

subroutine getElementArea(nPts,nDim,points,center,area,areaMag)
  use precision
  implicit none

  ! subroutine Variables
  integer(kind=intType)::nPts,nDim
  real(kind=realType),dimension(nPts,nDim)::points
  real(kind=realType),dimension(nDim)::center,area
  real(kind=realType)::areaMag

  ! Local variables
  integer(kind=intType):: i,j
  real(kind=realType),dimension(nPts,nDim):: radialVec
  real(kind=realType),dimension(nDim):: sumArea,cross
real(kind=realType),dimension(3):: v1, v2
  real(kind=realType):: tempArea

  ! begin execution
  ! compute the vector from the center to each point defining the element
  do i = 1,nPts
     radialVec(i,:) = points(i,:)-center(:)
  end do
  
  ! Now loop around element doing cross products to get directional area
  sumArea = 0
  do i = 1,nPts-1
     if (nDim .eq. 2) then
        ! Append a vector of zeros to the 2d pts
        v1 = 0.0
        v2 = 0.0
        v1(1:2) = radialVec(i,:)
        v1(1:2) = radialVec(i+1,:)
        ! now vectors exist in a hypothetical 3d plane
        call cross_product_3d(radialVec(i,:),radialVec(i+1,:),cross)
        sumArea = sumArea+ cross/2.0
     elseif( nDim .eq. 3) then
        call cross_product_3d(radialVec(i,:),radialVec(i+1,:),cross)
        sumArea = sumArea+ cross/2.0
     else
        print *,'Cross product not available for ',nDim,' dimensions...exiting'
        stop
     end if
  end do
  area = sumArea
 
  ! Get the magnitude of the area as well
  tempArea = 0
  do i = 1,nDim
     tempArea = tempArea+sumArea(i)*sumArea(i)
  end do
  areaMag = sqrt(tempArea)

end subroutine getElementArea
