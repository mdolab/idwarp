subroutine getSpherePts(x0, radius, vpts, fpts)
    ! Return the 20 points defining the nodes of a dodecahedron of
    ! radius 'radius' centered at point 'x0'. Also return the 12 points
    ! defining the center of each of the 12 faces.

    use constants
    implicit none

    real(kind=realType), intent(in) :: x0(3), radius
    real(kind=realType), intent(out) :: vpts(3, 20) ! vertex points
    real(kind=realType), intent(out) :: fpts(3, 12) ! face points
    integer(kind=intType) :: inds(5, 12), i, j
    real(kind=realType) :: ovrFive

    ! Nominal locations
    vpts(:, 1) = (/-one, -one, -one/)
    vpts(:, 2) = (/one, -one, -one/)
    vpts(:, 3) = (/-one, one, -one/)
    vpts(:, 4) = (/one, one, -one/)
    vpts(:, 5) = (/-one, -one, one/)
    vpts(:, 6) = (/one, -one, one/)
    vpts(:, 7) = (/-one, one, one/)
    vpts(:, 8) = (/one, one, one/)

    vpts(:, 9) = (/zero, -one / phi, -phi/)
    vpts(:, 10) = (/zero, one / phi, -phi/)
    vpts(:, 11) = (/zero, -one / phi, phi/)
    vpts(:, 12) = (/zero, one / phi, phi/)

    vpts(:, 13) = (/-one / phi, -phi, zero/)
    vpts(:, 14) = (/one / phi, -phi, zero/)
    vpts(:, 15) = (/-one / phi, phi, zero/)
    vpts(:, 16) = (/one / phi, phi, zero/)

    vpts(:, 17) = (/-phi, zero, -one / phi/)
    vpts(:, 18) = (/-phi, zero, one / phi/)
    vpts(:, 19) = (/phi, zero, -one / phi/)
    vpts(:, 20) = (/phi, zero, one / phi/)

    ! Now the indices
    inds(:, 1) = (/3, 4, 10, 15, 16/)
    inds(:, 2) = (/7, 8, 12, 15, 16/)
    inds(:, 3) = (/4, 8, 16, 19, 20/)
    inds(:, 4) = (/3, 7, 15, 17, 18/)

    inds(:, 5) = (/1, 3, 9, 10, 17/)
    inds(:, 6) = (/2, 4, 9, 10, 19/)
    inds(:, 7) = (/5, 7, 11, 12, 18/)
    inds(:, 8) = (/6, 8, 11, 12, 20/)

    inds(:, 9) = (/1, 2, 9, 13, 14/)
    inds(:, 10) = (/5, 6, 11, 13, 14/)
    inds(:, 11) = (/2, 6, 14, 19, 20/)
    inds(:, 12) = (/1, 5, 13, 17, 18/)

    ! Now scale by the desired radius...and the radius of nominal
    ! points is sqrt(3) and relocate
    do i = 1, 20
        vpts(:, i) = vpts(:, i) * radius / sqrt(three) + X0
    end do

    ! Now just average the nodes of the faces to get the face centers.
    fpts = zero
    do i = 1, 12
        do j = 1, 5
            fpts(:, i) = fpts(:, i) + vpts(:, inds(j, i))
        end do
    end do
    ovrFive = one / five
    fpts = fpts * ovrFive
end subroutine getSpherePts
