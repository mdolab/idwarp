! ====================================================================
! File: vectorUtils.f90
! Author: C.A.(Sandy) Mader
! Date Started: July 14, 2014
! Date Modified:

subroutine cross_product_3d(v1, v2, cross)
  use precision
  implicit none

  real(kind=realType), dimension(3), intent(in) :: v1, v2
  real(kind=realType), dimension(3), intent(out) :: cross
  
  cross(1) = (v1(2)*v2(3) - v1(3)*v2(2))
  cross(2) = (v1(3)*v2(1) - v1(1)*v2(3))
  cross(3) = (v1(1)*v2(2) - v1(2)*v2(1))

end subroutine cross_product_3d

subroutine getMag(V, mag)

  use constants
  implicit none

  ! Subroutine Variables
  real(kind=realType), dimension(3), intent(in) :: V
  real(kind=realType), intent(out) :: mag

  mag = sqrt(v(1)**2 + v(2)**2 + v(3)**2 + 1e-30)
  
end subroutine getMag

subroutine getRotationMatrix3d(v1, v2, Mi)
  use precision
  use constants
  implicit none

  ! Subroutine Variables
  real(kind=realType), dimension(3), intent(in) :: v1, v2
  real(kind=realType), dimension(3, 3), intent(out) :: Mi

  ! Local Variables
  real(kind=realType), dimension(3) :: axis
  real(kind=realType):: magV1, magV2, axisMag, angle
  real(kind=realType), dimension(3, 3):: axis_skewed
  logical :: flag
  call getMag(v1, magV1)
  call getMag(v2, magV2)
  
  ! Start by determining the rotation axis by getting the 
  ! cross product between v1, v2

  call cross_product_3d(v1,v2,axis)
  ! Now Normalize
  call getMag(axis,axisMag)
  flag = .False.
  if (axisMag <1.0e-8) then
     ! no rotation at this point, angle is 0
     angle = zero
     ! the axis doesn't matter so set to x
     axis = zero
     axis(1) = one
     flag = .True.
  else
     axis = axis/axisMag
     ! Now compute the rotation angle about that axis
     angle = acos(dot_product(v1/magv1,v2/magv2))
  end if
  ! Now that we have an axis and an angle,build the rotation Matrix

! A skew symmetric representation of the normalized axis 
  axis_skewed(1,1) = zero
  axis_skewed(1,2) = -axis(3)
  axis_skewed(1,3) = axis(2)
  axis_skewed(2,1) = axis(3)
  axis_skewed(2,2) = zero
  axis_skewed(2,3) = -axis(1)
  axis_skewed(3,1) = -axis(2)
  axis_skewed(3,2) = axis(1)
  axis_skewed(3,3) = zero

  ! identity matrix
  Mi = zero
  Mi(1, 1) = one
  Mi(2, 2) = one
  Mi(3, 3) = one
  
  ! Rodrigues formula for the rotation matrix 
  Mi = Mi + sin(angle)*axis_skewed + (one-cos(angle))*matmul(axis_skewed,axis_skewed)

end subroutine getRotationMatrix3d
