subroutine getElementProps(points, nPts, area, normal)
    use constants
    implicit none

    ! Subroutine Variables
    integer(kind=intType), intent(in) :: nPts
    real(kind=realType), dimension(3, nPts), intent(in) :: points
    real(kind=realType), dimension(3), intent(out) :: normal
    real(kind=realType), intent(out) :: area

    ! Local variables
    integer(kind=intType) :: i
    real(kind=realType), dimension(3, nPts) :: radialVec
    real(kind=realType), dimension(3) :: center, cross
    real(kind=realType), dimension(1) :: crossNorm

    ! Compute the element center:
    center(:) = zero
    do i = 1, nPts
        center(:) = center(:) + points(:, i)
    end do
    center = center / nPts

    ! compute the vector from the center to each point defining the element
    do i = 1, nPts
        radialVec(:, i) = points(:, i) - center
    end do

    ! Now loop around element doing cross products to get directional area
    normal = zero
    area = zero
    do i = 1, nPts
        if (i < nPts) then
            call cross_product_3d(radialVec(:, i), radialVec(:, i + 1), cross)
        else
            call cross_product_3d(radialVec(:, i), radialVec(:, 1), cross)
        end if
        crossNorm(1) = sqrt(cross(1)**2 + cross(2)**2 + cross(3)**2 + 1e-15)
        area = area + half * crossNorm(1)
        normal = normal + cross / crossNorm(1)
    end do

    ! Also make the normal a unit lenth
    normal = normal / sqrt(normal(1)**2 + normal(2)**2 + normal(3)**2)

end subroutine getElementProps
