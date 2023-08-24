   !        Generated by TAPENADE     (INRIA, Ecuador team)
   !  Tapenade 3.16 (develop) - 11 Aug 2023 14:53
   !
   !  Differentiation of getrotationmatrix3d in reverse (adjoint) mode (with options noISIZE i4 dr8 r8):
   !   gradient     of useful results: v2 mi
   !   with respect to varying inputs: v2
   SUBROUTINE GETROTATIONMATRIX3D_B(v1, v2, v2b, mi, mib)
   USE PRECISION
   USE CONSTANTS
   IMPLICIT NONE
   ! Subroutine Variables
   REAL(kind=realtype), DIMENSION(3), INTENT(IN) :: v1, v2
   REAL(kind=realtype), DIMENSION(3) :: v2b
   REAL(kind=realtype), DIMENSION(3, 3) :: mi
   REAL(kind=realtype), DIMENSION(3, 3) :: mib
   ! Local Variables
   REAL(kind=realtype), DIMENSION(3) :: axis, vv1, vv2
   REAL(kind=realtype), DIMENSION(3) :: axisb, vv2b
   REAL(kind=realtype) :: magv1, magv2, axismag, angle, arg
   REAL(kind=realtype) :: magv2b, axismagb, angleb, argb
   REAL(kind=realtype), DIMENSION(3, 3) :: a, c
   REAL(kind=realtype), DIMENSION(3, 3) :: ab, cb
   REAL(kind=realtype), PARAMETER :: tol=1.4901161193847656e-08
   INTRINSIC MIN
   INTRINSIC ACOS
   INTRINSIC SIN
   INTRINSIC COS
   REAL(kind=realtype), DIMENSION(3) :: v1b
   INTEGER :: branch
   CALL GETMAG(v1, magv1)
   CALL GETMAG(v2, magv2)
   ! Start by determining the rotation axis by getting the
   ! cross product between v1, v2
   CALL CROSS_PRODUCT_3D(v1, v2, axis)
   ! Now Normalize
   CALL GETMAG(axis, axismag)
   ! When axisMag is less that sqrt(eps), the acos 'arg' value will be
   ! exactly one which will give a nan in complex mode.
   IF (axismag .LT. tol) THEN
   ! no rotation at this point, angle is 0
   angle = zero
   ! the axis doesn't matter so set to x
   CALL PUSHREAL8ARRAY(axis, realtype*3/8)
   axis = zero
   axis(1) = one
   CALL PUSHCONTROL1B(0)
   ELSE
   CALL PUSHREAL8ARRAY(axis, realtype*3/8)
   axis = axis/axismag
   ! Now compute the rotation angle about that axis
   vv1 = v1/magv1
   vv2 = v2/magv2
   IF (one .GT. vv1(1)*vv2(1) + vv1(2)*vv2(2) + vv1(3)*vv2(3)) THEN
   arg = vv1(1)*vv2(1) + vv1(2)*vv2(2) + vv1(3)*vv2(3)
   CALL PUSHCONTROL1B(0)
   ELSE
   arg = one
   CALL PUSHCONTROL1B(1)
   END IF
   angle = ACOS(arg)
   CALL PUSHCONTROL1B(1)
   END IF
   ! Now that we have an axis and an angle,build the rotation Matrix
   ! A skew symmetric representation of the normalized axis
   a(1, 1) = zero
   a(1, 2) = -axis(3)
   a(1, 3) = axis(2)
   a(2, 1) = axis(3)
   a(2, 2) = zero
   a(2, 3) = -axis(1)
   a(3, 1) = -axis(2)
   a(3, 2) = axis(1)
   a(3, 3) = zero
   !C = A*A
   c(1, 1) = a(1, 1)*a(1, 1) + a(1, 2)*a(2, 1) + a(1, 3)*a(3, 1)
   c(1, 2) = a(1, 1)*a(1, 2) + a(1, 2)*a(2, 2) + a(1, 3)*a(3, 2)
   c(1, 3) = a(1, 1)*a(1, 3) + a(1, 2)*a(2, 3) + a(1, 3)*a(3, 3)
   c(2, 1) = a(2, 1)*a(1, 1) + a(2, 2)*a(2, 1) + a(2, 3)*a(3, 1)
   c(2, 2) = a(2, 1)*a(1, 2) + a(2, 2)*a(2, 2) + a(2, 3)*a(3, 2)
   c(2, 3) = a(2, 1)*a(1, 3) + a(2, 2)*a(2, 3) + a(2, 3)*a(3, 3)
   c(3, 1) = a(3, 1)*a(1, 1) + a(3, 2)*a(2, 1) + a(3, 3)*a(3, 1)
   c(3, 2) = a(3, 1)*a(1, 2) + a(3, 2)*a(2, 2) + a(3, 3)*a(3, 2)
   c(3, 3) = a(3, 1)*a(1, 3) + a(3, 2)*a(2, 3) + a(3, 3)*a(3, 3)
   ! Rodrigues formula for the rotation matrix
   ab = 0.0_8
   cb = 0.0_8
   angleb = COS(angle)*SUM(a*mib) + SIN(angle)*SUM(c*mib)
   ab = SIN(angle)*mib
   cb = (one-COS(angle))*mib
   ab(3, 1) = ab(3, 1) + a(1, 3)*cb(3, 3) + a(1, 2)*cb(3, 2) + (a(1, 1)+a&
   &   (3, 3))*cb(3, 1) + a(2, 3)*cb(2, 1) + a(1, 3)*cb(1, 1)
   ab(1, 3) = ab(1, 3) + a(3, 1)*cb(3, 3) + a(2, 1)*cb(2, 3) + (a(1, 1)+a&
   &   (3, 3))*cb(1, 3) + a(3, 2)*cb(1, 2) + a(3, 1)*cb(1, 1)
   ab(3, 2) = ab(3, 2) + a(2, 3)*cb(3, 3) + (a(2, 2)+a(3, 3))*cb(3, 2) + &
   &   a(2, 1)*cb(3, 1) + a(2, 3)*cb(2, 2) + a(1, 3)*cb(1, 2)
   ab(2, 3) = ab(2, 3) + a(3, 2)*cb(3, 3) + (a(2, 2)+a(3, 3))*cb(2, 3) + &
   &   a(3, 2)*cb(2, 2) + a(3, 1)*cb(2, 1) + a(1, 2)*cb(1, 3)
   ab(1, 2) = ab(1, 2) + a(3, 1)*cb(3, 2) + a(2, 1)*cb(2, 2) + a(2, 3)*cb&
   &   (1, 3) + (a(1, 1)+a(2, 2))*cb(1, 2) + a(2, 1)*cb(1, 1)
   ab(1, 1) = ab(1, 1) + a(3, 1)*cb(3, 1) + a(2, 1)*cb(2, 1) + a(1, 3)*cb&
   &   (1, 3) + a(1, 2)*cb(1, 2) + 2*a(1, 1)*cb(1, 1)
   ab(2, 1) = ab(2, 1) + a(3, 2)*cb(3, 1) + a(1, 3)*cb(2, 3) + a(1, 2)*cb&
   &   (2, 2) + (a(1, 1)+a(2, 2))*cb(2, 1) + a(1, 2)*cb(1, 1)
   ab(3, 3) = 0.0_8
   cb(3, 3) = 0.0_8
   cb(3, 1) = 0.0_8
   cb(1, 3) = 0.0_8
   axisb = 0.0_8
   axisb(1) = axisb(1) + ab(3, 2) - ab(2, 3)
   ab(3, 2) = 0.0_8
   axisb(2) = axisb(2) + ab(1, 3) - ab(3, 1)
   ab(3, 1) = 0.0_8
   ab(2, 3) = 0.0_8
   ab(2, 2) = 0.0_8
   cb(3, 2) = 0.0_8
   cb(2, 3) = 0.0_8
   cb(2, 2) = 0.0_8
   cb(2, 1) = 0.0_8
   cb(1, 2) = 0.0_8
   axisb(3) = axisb(3) + ab(2, 1) - ab(1, 2)
   ab(2, 1) = 0.0_8
   ab(1, 3) = 0.0_8
   CALL POPCONTROL1B(branch)
   IF (branch .EQ. 0) THEN
   CALL POPREAL8ARRAY(axis, realtype*3/8)
   magv2b = 0.0_8
   axisb = 0.0_8
   axismagb = 0.0_8
   ELSE
   IF (arg .EQ. 1.0 .OR. arg .EQ. (-1.0)) THEN
   argb = 0.0_8
   ELSE
   argb = -(angleb/SQRT(1.0-arg**2))
   END IF
   CALL POPCONTROL1B(branch)
   IF (branch .EQ. 0) THEN
   vv2b = 0.0_8
   vv2b(1) = vv2b(1) + vv1(1)*argb
   vv2b(2) = vv2b(2) + vv1(2)*argb
   vv2b(3) = vv2b(3) + vv1(3)*argb
   ELSE
   vv2b = 0.0_8
   END IF
   v2b = v2b + vv2b/magv2
   magv2b = -(SUM(v2*vv2b)/magv2**2)
   CALL POPREAL8ARRAY(axis, realtype*3/8)
   axismagb = -(SUM(axis*axisb)/axismag**2)
   axisb = axisb/axismag
   END IF
   CALL GETMAG_B(axis, axisb, axismag, axismagb)
   v1b = 0.0_8
   CALL CROSS_PRODUCT_3D_B(v1, v1b, v2, v2b, axis, axisb)
   CALL GETMAG_B(v2, v2b, magv2, magv2b)
   END SUBROUTINE GETROTATIONMATRIX3D_B
      !  Differentiation of cross_product_3d in reverse (adjoint) mode (with options noISIZE i4 dr8 r8):
   !   gradient     of useful results: cross v1 v2
   !   with respect to varying inputs: cross v1 v2
   ! ====================================================================
   ! File: vectorUtils.f90
   ! Author: C.A.(Sandy) Mader
   ! Date Started: July 14, 2014
   ! Date Modified:
   SUBROUTINE CROSS_PRODUCT_3D_B(v1, v1b, v2, v2b, cross, crossb)
   USE CONSTANTS
   IMPLICIT NONE
   REAL(kind=realtype), DIMENSION(3), INTENT(IN) :: v1, v2
   REAL(kind=realtype), DIMENSION(3) :: v1b, v2b
   REAL(kind=realtype), DIMENSION(3) :: cross
   REAL(kind=realtype), DIMENSION(3) :: crossb
   v1b(1) = v1b(1) + v2(2)*crossb(3) - v2(3)*crossb(2)
   v2b(2) = v2b(2) + v1(1)*crossb(3) - v1(3)*crossb(1)
   v1b(2) = v1b(2) + v2(3)*crossb(1) - v2(1)*crossb(3)
   v2b(1) = v2b(1) + v1(3)*crossb(2) - v1(2)*crossb(3)
   crossb(3) = 0.0_8
   v1b(3) = v1b(3) + v2(1)*crossb(2) - v2(2)*crossb(1)
   v2b(3) = v2b(3) + v1(2)*crossb(1) - v1(1)*crossb(2)
   crossb(2) = 0.0_8
   crossb(1) = 0.0_8
   END SUBROUTINE CROSS_PRODUCT_3D_B
      !  Differentiation of getmag in reverse (adjoint) mode (with options noISIZE i4 dr8 r8):
   !   gradient     of useful results: v mag
   !   with respect to varying inputs: v
   SUBROUTINE GETMAG_B(v, vb, mag, magb)
   USE CONSTANTS
   IMPLICIT NONE
   ! Subroutine Variables
   REAL(kind=realtype), DIMENSION(3), INTENT(IN) :: v
   REAL(kind=realtype), DIMENSION(3) :: vb
   REAL(kind=realtype) :: mag
   REAL(kind=realtype) :: magb
   INTRINSIC SQRT
   REAL(kind=realtype) :: tempb
   IF (v(1)**2 + v(2)**2 + v(3)**2 + 1e-30 .EQ. 0.0_8) THEN
   tempb = 0.0_8
   ELSE
   tempb = magb/(2.0*SQRT(v(1)**2+v(2)**2+v(3)**2+1e-30))
   END IF
   vb(1) = vb(1) + 2*v(1)*tempb
   vb(2) = vb(2) + 2*v(2)*tempb
   vb(3) = vb(3) + 2*v(3)*tempb
   END SUBROUTINE GETMAG_B
   