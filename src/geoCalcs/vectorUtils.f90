! ====================================================================
! File: vectorUtils.f90
! Author: C.A.(Sandy) Mader
! Date Started: July 14, 2014
! Date Modified:


subroutine cross_product_3d(v1,v2,cross)
  use precision
  implicit none

  real(kind=realType),dimension(3)::v1,v2,cross
  
  cross(1) = half*(v1(2)*v2(3) - v1(3)*v2(2))
  cross(2) = half*(v1(3)*v2(1) - v1(1)*v2(3))
  cross(3) = half*(v1(1)*v2(2) - v1(2)*v2(1))
end subroutine cross_product_3d
