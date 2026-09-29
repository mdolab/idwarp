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

!> @brief Compute the rotation matrix that maps one vector to another
!!
!! If the two vectors are not unit vectors then the rotation matrix is actually for the transformation between unit
!! vectors in the directions of the two input vectors
!!
!! @param[in]  v1 Source vector
!! @param[in]  v2 Target vector
!! @param[out] R Rotation matrix
subroutine getRotationMatrix3d(v1, v2, R)
    use precision
    use constants
    implicit none

    ! Subroutine Variables
    real(kind=realType), dimension(3), intent(in) :: v1, v2
    real(kind=realType), dimension(3, 3), intent(out) :: R

    ! Local Variables
    real(kind=realType), dimension(3) :: v1Unit, v2Unit, v
    real(kind=realType) :: v1Mag, v2Mag, c, denom
    real(kind=realType), dimension(3, 3) :: skew, skew2
    real(kind=realType), parameter :: dtol = 1.0e-10

    ! Based on Rodrigues formula, the rotation matrix between two unit vectors is:
    ! R = I + skew(v) + 1 / (1 + c) * skew(v)^2
    ! Where:
    ! - v = cross(v1, v2)
    ! - skew(v) is the cross product matrix of v, such that skew(v)*x = cross(v,x)
    ! - c = dot(v1, v2)
    ! Source: https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d/476311#476311

    ! We use dtol to clip the value of 1/(1+c) when v1 and v2 are close to anti-parallel since in this case c=-1 and
    ! 1 / (1 + c) is undefined. This reflects the fact that there is no unique rotation axis for the anti-parallel case

    ! Normalize input vectors first so the above formula holds
    call getMag(v1, v1Mag)
    v1Unit = v1 / v1Mag
    call getMag(v2, v2Mag)
    v2Unit = v2 / v2Mag
    call cross_product_3d(v1Unit, v2Unit, v)
    c = v1Unit(1) * v2Unit(1) + v1Unit(2) * v2Unit(2) + v1Unit(3) * v2Unit(3)
    denom = one + c
    if (denom < dtol) then
        denom = dtol
    end if

    ! skew(v)
    skew(1, 1) = zero
    skew(1, 2) = -v(3)
    skew(1, 3) = v(2)
    skew(2, 1) = v(3)
    skew(2, 2) = zero
    skew(2, 3) = -v(1)
    skew(3, 1) = -v(2)
    skew(3, 2) = v(1)
    skew(3, 3) = zero

    ! skew(v)^2 = v v^T - (v.v) I
    skew2(1, 1) = -v(2) * v(2) - v(3) * v(3)
    skew2(1, 2) = v(1) * v(2)
    skew2(1, 3) = v(1) * v(3)
    skew2(2, 1) = v(1) * v(2)
    skew2(2, 2) = -v(1) * v(1) - v(3) * v(3)
    skew2(2, 3) = v(2) * v(3)
    skew2(3, 1) = v(1) * v(3)
    skew2(3, 2) = v(2) * v(3)
    skew2(3, 3) = -v(1) * v(1) - v(2) * v(2)

    R = zero
    R(1, 1) = one
    R(2, 2) = one
    R(3, 3) = one

    R = R + skew + skew2 / denom

end subroutine getRotationMatrix3d
