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
    real(kind=realType), dimension(3) :: a
    real(kind=realType) :: magV1, magV2, s, dot, denom2
    real(kind=realType), dimension(3, 3) :: P, C
    real(kind=realType), parameter :: dtol = 1.0e-10

    ! Rotation taking the direction of v1 to the direction of v2, in a
    ! smooth branch-free form.  With theta the angle between v1 and v2,
    ! a = v1 x v2 and s = |v1||v2|, the two axis-angle pieces of the
    ! Rodrigues formula reduce EXACTLY (for 0 <= theta < pi) to
    !     sin(theta) * Ahat        = skew(a) / s
    !     (1-cos(theta)) * Ahat^2  = skew(a)^2 / (s^2 + s*(v1.v2))
    ! so no normalized axis, no acos and no angle are ever formed.  The
    ! axis-angle version this replaces needed an axisMag < tol branch at
    ! theta = 0, and BOTH AD derivatives (_b and _d) of that branch were
    ! identically zero while the true derivative there is finite: every
    ! surface node of an UNDEFORMED mesh evaluates at exactly theta = 0,
    ! so warpDeriv/warpDerivFwd silently dropped the entire rotation
    ! contribution when linearized at the baseline configuration.  This
    ! form has a removable limit at theta = 0 (a = 0 makes the rotation
    ! terms vanish smoothly) and differentiates correctly there by
    ! construction.  Only theta -> pi (anti-parallel normals, a folded
    ! surface) is guarded, where the rotation genuinely is not unique.

    call getMag(v1, magV1)
    call getMag(v2, magV2)
    call cross_product_3d(v1, v2, a)
    s = magV1 * magV2
    dot = v1(1) * v2(1) + v1(2) * v2(2) + v1(3) * v2(3)
    denom2 = s * s + s * dot
    if (denom2 < dtol * s * s) then
        denom2 = dtol * s * s
    end if

    ! P = skew(a)
    P(1, 1) = zero
    P(1, 2) = -a(3)
    P(1, 3) = a(2)
    P(2, 1) = a(3)
    P(2, 2) = zero
    P(2, 3) = -a(1)
    P(3, 1) = -a(2)
    P(3, 2) = a(1)
    P(3, 3) = zero

    ! C = P*P = a a^T - (a.a) I
    C(1, 1) = -a(2) * a(2) - a(3) * a(3)
    C(1, 2) = a(1) * a(2)
    C(1, 3) = a(1) * a(3)
    C(2, 1) = a(1) * a(2)
    C(2, 2) = -a(1) * a(1) - a(3) * a(3)
    C(2, 3) = a(2) * a(3)
    C(3, 1) = a(1) * a(3)
    C(3, 2) = a(2) * a(3)
    C(3, 3) = -a(1) * a(1) - a(2) * a(2)

    Mi = zero
    Mi(1, 1) = one
    Mi(2, 2) = one
    Mi(3, 3) = one

    Mi = Mi + P / s + C / denom2

end subroutine getRotationMatrix3d
