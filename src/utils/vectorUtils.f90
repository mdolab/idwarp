! ====================================================================
! File: vectorUtils.f90
! Author: C.A.(Sandy) Mader
! Date Started: July 14, 2014
! Date Modified:

subroutine cross_product_3d(v1, v2, cross)
    use constants
    implicit none

    real(kind=realType), dimension(3), intent(in) :: v1, v2
    real(kind=realType), dimension(3), intent(out) :: cross
    cross(1) = (v1(2) * v2(3) - v1(3) * v2(2))
    cross(2) = (v1(3) * v2(1) - v1(1) * v2(3))
    cross(3) = (v1(1) * v2(2) - v1(2) * v2(1))

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
    real(kind=realType), dimension(3) :: axis, vv1, vv2
    real(kind=realType) :: magV1, magV2, axisMag, angle, arg
    real(kind=realType), dimension(3, 3) :: A, C
    real(kind=realType), parameter :: tol = 1.4901161193847656e-08

    call getMag(v1, magV1)
    call getMag(v2, magV2)

    ! Start by determining the rotation axis by getting the
    ! cross product between v1, v2

    call cross_product_3d(v1, v2, axis)
    ! Now Normalize
    call getMag(axis, axisMag)

    ! When axisMag is less that sqrt(eps), the acos 'arg' value will be
    ! exactly one which will give a nan in complex mode.
    if (axisMag < tol) then
        ! no rotation at this point, angle is 0
        angle = zero
        ! the axis doesn't matter so set to x
        axis = zero
        axis(1) = one
    else
        axis = axis / axisMag
        ! Now compute the rotation angle about that axis
        vv1 = v1 / magv1
        vv2 = v2 / magv2
        arg = min(one, vv1(1) * vv2(1) + vv1(2) * vv2(2) + vv1(3) * vv2(3))
        angle = acos(arg)
    end if
    ! Now that we have an axis and an angle,build the rotation Matrix
    ! A skew symmetric representation of the normalized axis
    A(1, 1) = zero
    A(1, 2) = -axis(3)
    A(1, 3) = axis(2)
    A(2, 1) = axis(3)
    A(2, 2) = zero
    A(2, 3) = -axis(1)
    A(3, 1) = -axis(2)
    A(3, 2) = axis(1)
    A(3, 3) = zero

    !C = A*A
    C(1, 1) = A(1, 1) * A(1, 1) + A(1, 2) * A(2, 1) + A(1, 3) * A(3, 1)
    C(1, 2) = A(1, 1) * A(1, 2) + A(1, 2) * A(2, 2) + A(1, 3) * A(3, 2)
    C(1, 3) = A(1, 1) * A(1, 3) + A(1, 2) * A(2, 3) + A(1, 3) * A(3, 3)
    C(2, 1) = A(2, 1) * A(1, 1) + A(2, 2) * A(2, 1) + A(2, 3) * A(3, 1)
    C(2, 2) = A(2, 1) * A(1, 2) + A(2, 2) * A(2, 2) + A(2, 3) * A(3, 2)
    C(2, 3) = A(2, 1) * A(1, 3) + A(2, 2) * A(2, 3) + A(2, 3) * A(3, 3)
    C(3, 1) = A(3, 1) * A(1, 1) + A(3, 2) * A(2, 1) + A(3, 3) * A(3, 1)
    C(3, 2) = A(3, 1) * A(1, 2) + A(3, 2) * A(2, 2) + A(3, 3) * A(3, 2)
    C(3, 3) = A(3, 1) * A(1, 3) + A(3, 2) * A(2, 3) + A(3, 3) * A(3, 3)

    ! Rodrigues formula for the rotation matrix
    Mi = zero
    Mi(1, 1) = one
    Mi(2, 2) = one
    Mi(3, 3) = one

    Mi = Mi + sin(angle) * A + (one - cos(angle)) * C

end subroutine getRotationMatrix3d
