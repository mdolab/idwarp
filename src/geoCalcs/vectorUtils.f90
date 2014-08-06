! ====================================================================
! File: vectorUtils.f90
! Author: C.A.(Sandy) Mader
! Date Started: July 14, 2014
! Date Modified:


subroutine cross_product_3d(v1,v2,cross)
  use precision
  use constants
  implicit none

  real(kind=realType),dimension(3)::v1,v2,cross
  
  cross(1) = half*(v1(2)*v2(3) - v1(3)*v2(2))
  cross(2) = half*(v1(3)*v2(1) - v1(1)*v2(3))
  cross(3) = half*(v1(1)*v2(2) - v1(2)*v2(1))
end subroutine cross_product_3d

! subroutine getDistance(nDim,r,ri,dist)
!   use precision
!   use constants
!   implicit none

!   ! Subroutine Variables
!   integer(kind=intType)::nDim
!   real(kind=realType),dimension(nDim)::r,ri
!   real(kind=realType):: dist

!   ! Local Variables
!   real(kind=realType),dimension(nDim)::dx

!   dx = r-ri
!   call getMag(nDim,dx,dist)

! end subroutine getDistance

subroutine getDistance(r,ri,dist)
  use precision
  use constants
  implicit none

  ! Subroutine Variables
  integer(kind=intType)::nDim
  real(kind=realType),dimension(3)::r,ri
  real(kind=realType):: dist

  ! Local Variables
  real(kind=realType),dimension(3)::dx

  dx = r-ri
  call getMag(dx,dist)

end subroutine getDistance


! subroutine getMag(nDim,V,mag)
!   use precision
!   use constants
!   implicit none

!   ! Subroutine Variables
!   integer(kind=intType)::nDim
!   real(kind=realType),dimension(nDim)::V
!   real(kind=realType):: mag

!   ! Local Variables
!   integer(kind=intType)::i
!   real(kind=realType)::sum

!   sum=0
!   do i = 1,nDim
!      sum = sum+V(i)*V(i)
!   end do
!   mag = sqrt(sum)
  
! end subroutine getMag
subroutine getMag(V,mag)
  use precision
  use constants
  implicit none

  ! Subroutine Variables
  real(kind=realType),dimension(3)::V
  real(kind=realType):: mag

  ! Local Variables
  ! integer(kind=intType)::i
  ! real(kind=realType)::sum

  ! sum=0
  ! do i = 1,3
  !    sum = sum+V(i)*V(i)
  ! end do
  ! mag = sqrt(sum)
  mag = sqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))
  
end subroutine getMag

subroutine getRotationMatrix3d(v1,v2,Mi)
  use precision
  use constants
  implicit none

  ! Subroutine Variables
  real(kind=realType),dimension(3)::v1,v2
  real(kind=realType),dimension(3,3):: Mi

  ! Local Variables
  real(kind=realType),dimension(3)::axis
  real(kind=realType)::magV1,magV2,axisMag,angle
  real(kind=realType),dimension(3,3)::eye, axis_skewed
  integer(kind=intType)::i

  ! Begin Execution
  !print *,'v1,v2',v1,v2
  call getMag(v1,magV1)
  call getMag(v2,magV2)
  
  !print *,'mags',magv1,magv2

  ! Start by determining the rotation axis by getting the 
  ! cross product between v1, v2

  call cross_product_3d(v1,v2,axis)
  ! Now Normalize
  call getMag(axis,axisMag)
  !print *,'axis',axis,'mag', axisMag
  if (axisMag <1.0e-15) then
     ! no rotation at this point, angle is 0
     angle = 0
     ! the axis doesn't matter so set to x
     axis = 0
     axis(1) = 1.0
  else
     axis = axis/axisMag

     ! Now compute the rotation angle about that axis
     angle = acos(dot_product(v1/magv1,v2/magv2))
     !print *,'angle',angle,'ax',axis
  end if
  ! Now that we have an axis and an angle,build the rotation Matrix

  ! A skew symmetric representation of the normalized axis 
  axis_skewed(1,1) = 0.
  axis_skewed(1,2) = -axis(3)
  axis_skewed(1,3) = axis(2)
  axis_skewed(2,1) = axis(3)
  axis_skewed(2,2) = 0.
  axis_skewed(2,3) = -axis(1)
  axis_skewed(3,1) = -axis(2)
  axis_skewed(3,2) = axis(1)
  axis_skewed(3,3) = 0.

  ! identity matrix
  eye = 0
  do i = 1,3
     eye(i,i) = 1.0
  end do

  ! Rodrigues formula for the rotation matrix 
  
  Mi = eye + sin(angle)*axis_skewed + (1-cos(angle))*matmul(axis_skewed,axis_skewed)
  !print *,'Rot:Mi',Mi,'a',angle,'ax',axis
  ! stop
end subroutine getRotationMatrix3d
