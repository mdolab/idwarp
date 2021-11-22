   !        Generated by TAPENADE     (INRIA, Ecuador team)
   !  Tapenade 3.16 (master) -  9 Oct 2020 17:47
   !
   !  Differentiation of getelementprops in forward (tangent) mode (with options noISIZE i4 dr8 r8):
   !   variations   of useful results: area normal
   !   with respect to varying inputs: points
   SUBROUTINE GETELEMENTPROPS_D(points, pointsb, npts, area, areab, normal&
   & , normalb)
   USE CONSTANTS
   IMPLICIT NONE
   ! Subroutine Variables
   INTEGER(kind=inttype), INTENT(IN) :: npts
   REAL(kind=realtype), DIMENSION(3, npts), INTENT(IN) :: points
   REAL(kind=realtype), DIMENSION(3, npts), INTENT(IN) :: pointsb
   REAL(kind=realtype), DIMENSION(3), INTENT(OUT) :: normal
   REAL(kind=realtype), DIMENSION(3), INTENT(OUT) :: normalb
   REAL(kind=realtype), INTENT(OUT) :: area
   REAL(kind=realtype), INTENT(OUT) :: areab
   ! Local variables
   INTEGER(kind=inttype) :: i
   REAL(kind=realtype), DIMENSION(3, npts) :: radialvec
   REAL(kind=realtype), DIMENSION(3, npts) :: radialvecb
   REAL(kind=realtype), DIMENSION(3) :: center, cross
   REAL(kind=realtype), DIMENSION(3) :: centerb, crossb
   REAL(kind=realtype), DIMENSION(1) :: crossnorm
   REAL(kind=realtype), DIMENSION(1) :: crossnormb
   INTRINSIC SQRT
   REAL(kind=realtype) :: arg1
   REAL(kind=realtype) :: arg1b
   REAL(kind=realtype) :: result1
   REAL(kind=realtype) :: result1b
   REAL(kind=realtype) :: temp
   ! Compute the element center:
   center(:) = zero
   centerb = 0.0_8
   DO i=1,npts
   centerb(:) = centerb(:) + pointsb(:, i)
   center(:) = center(:) + points(:, i)
   END DO
   centerb = centerb/npts
   center = center/npts
   radialvecb = 0.0_8
   ! compute the vector from the center to each point defining the element
   DO i=1,npts
   radialvecb(:, i) = pointsb(:, i) - centerb
   radialvec(:, i) = points(:, i) - center
   END DO
   ! Now loop around element doing cross products to get directional area
   normal = zero
   area = zero
   areab = 0.0_8
   normalb = 0.0_8
   crossb = 0.0_8
   crossnormb = 0.0_8
   DO i=1,npts
   IF (i .LT. npts) THEN
   CALL CROSS_PRODUCT_3D_D(radialvec(:, i), radialvecb(:, i), &
   &                       radialvec(:, i+1), radialvecb(:, i+1), cross, &
   &                       crossb)
   ELSE
   CALL CROSS_PRODUCT_3D_D(radialvec(:, i), radialvecb(:, i), &
   &                       radialvec(:, 1), radialvecb(:, 1), cross, crossb&
   &                      )
   END IF
   arg1b = 2*cross(1)*crossb(1) + 2*cross(2)*crossb(2) + 2*cross(3)*&
   &     crossb(3)
   arg1 = cross(1)**2 + cross(2)**2 + cross(3)**2 + 1e-15
   temp = SQRT(arg1)
   IF (arg1 .EQ. 0.0_8) THEN
   crossnormb(1) = 0.0_8
   ELSE
   crossnormb(1) = arg1b/(2.0*temp)
   END IF
   crossnorm(1) = temp
   areab = areab + half*crossnormb(1)
   area = area + half*crossnorm(1)
   normalb = normalb + (crossb-cross*crossnormb(1)/crossnorm(1))/&
   &     crossnorm(1)
   normal = normal + cross/crossnorm(1)
   END DO
   ! Also make the normal a unit lenth
   arg1b = 2*normal(1)*normalb(1) + 2*normal(2)*normalb(2) + 2*normal(3)*&
   &   normalb(3)
   arg1 = normal(1)**2 + normal(2)**2 + normal(3)**2
   temp = SQRT(arg1)
   IF (arg1 .EQ. 0.0_8) THEN
   result1b = 0.0_8
   ELSE
   result1b = arg1b/(2.0*temp)
   END IF
   result1 = temp
   normalb = (normalb-normal*result1b/result1)/result1
   normal = normal/result1
   END SUBROUTINE GETELEMENTPROPS_D
   