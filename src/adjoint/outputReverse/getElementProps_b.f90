   !        Generated by TAPENADE     (INRIA, Ecuador team)
   !  Tapenade 3.16 (master) -  9 Oct 2020 17:47
   !
   !  Differentiation of getelementprops in reverse (adjoint) mode (with options noISIZE i4 dr8 r8):
   !   gradient     of useful results: area points normal
   !   with respect to varying inputs: points
   SUBROUTINE GETELEMENTPROPS_B(points, pointsb, npts, area, areab, normal&
   & , normalb)
   USE CONSTANTS
   IMPLICIT NONE
   ! Subroutine Variables
   INTEGER(kind=inttype), INTENT(IN) :: npts
   REAL(kind=realtype), DIMENSION(3, npts), INTENT(IN) :: points
   REAL(kind=realtype), DIMENSION(3, npts) :: pointsb
   REAL(kind=realtype), DIMENSION(3) :: normal
   REAL(kind=realtype), DIMENSION(3) :: normalb
   REAL(kind=realtype) :: area
   REAL(kind=realtype) :: areab
   ! Local variables
   INTEGER(kind=inttype) :: i
   REAL(kind=realtype), DIMENSION(3, npts) :: radialvec
   REAL(kind=realtype), DIMENSION(3, npts) :: radialvecb
   REAL(kind=realtype), DIMENSION(3) :: center, cross
   REAL(kind=realtype), DIMENSION(3) :: centerb, crossb
   REAL(kind=realtype), DIMENSION(1) :: crossnorm
   REAL(kind=realtype), DIMENSION(1) :: crossnormb
   INTRINSIC SQRT
   REAL(kind=realtype) :: temp
   REAL(kind=realtype) :: tempb
   REAL(kind=realtype), DIMENSION(3) :: tmp
   REAL(kind=realtype), DIMENSION(3) :: tmpb
   REAL(kind=realtype) :: temp0
   INTEGER :: branch
   ! Compute the element center:
   center(:) = zero
   DO i=1,npts
   center(:) = center(:) + points(:, i)
   END DO
   center = center/npts
   ! compute the vector from the center to each point defining the element
   DO i=1,npts
   radialvec(:, i) = points(:, i) - center
   END DO
   ! Now loop around element doing cross products to get directional area
   CALL PUSHREAL8ARRAY(normal, realtype*3/8)
   normal = zero
   DO i=1,npts
   IF (i .LT. npts) THEN
   CALL PUSHREAL8ARRAY(cross, realtype*3/8)
   CALL CROSS_PRODUCT_3D(radialvec(:, i), radialvec(:, i+1), cross)
   CALL PUSHCONTROL1B(0)
   ELSE
   CALL PUSHREAL8ARRAY(cross, realtype*3/8)
   CALL CROSS_PRODUCT_3D(radialvec(:, i), radialvec(:, 1), cross)
   CALL PUSHCONTROL1B(1)
   END IF
   CALL PUSHREAL8ARRAY(crossnorm(1), realtype/8)
   crossnorm(1) = SQRT(cross(1)**2 + cross(2)**2 + cross(3)**2 + 1e-15)
   normal = normal + cross/crossnorm(1)
   END DO
   temp = normal(1)*normal(1) + normal(2)*normal(2) + normal(3)*normal(3)
   temp0 = SQRT(temp)
   tmpb = normalb
   normalb = tmpb/temp0
   IF (temp .EQ. 0.0_8) THEN
   tempb = 0.0_8
   ELSE
   tempb = -(SUM(normal*tmpb)/(2.0*temp0**3))
   END IF
   normalb(1) = normalb(1) + 2*normal(1)*tempb
   normalb(2) = normalb(2) + 2*normal(2)*tempb
   normalb(3) = normalb(3) + 2*normal(3)*tempb
   crossb = 0.0_8
   radialvecb = 0.0_8
   crossnormb = 0.0_8
   DO i=npts,1,-1
   crossb = crossb + normalb/crossnorm(1)
   crossnormb(1) = crossnormb(1) + half*areab - SUM(cross*normalb)/&
   &     crossnorm(1)**2
   CALL POPREAL8ARRAY(crossnorm(1), realtype/8)
   IF (cross(1)**2 + cross(2)**2 + cross(3)**2 + 1e-15 .EQ. 0.0_8) THEN
   tempb = 0.0_8
   ELSE
   tempb = crossnormb(1)/(2.0*SQRT(cross(1)**2+cross(2)**2+cross(3)**&
   &       2+1e-15))
   END IF
   crossnormb(1) = 0.0_8
   crossb(1) = crossb(1) + 2*cross(1)*tempb
   crossb(2) = crossb(2) + 2*cross(2)*tempb
   crossb(3) = crossb(3) + 2*cross(3)*tempb
   CALL POPCONTROL1B(branch)
   IF (branch .EQ. 0) THEN
   CALL POPREAL8ARRAY(cross, realtype*3/8)
   CALL CROSS_PRODUCT_3D_B(radialvec(:, i), radialvecb(:, i), &
   &                       radialvec(:, i+1), radialvecb(:, i+1), cross, &
   &                       crossb)
   ELSE
   CALL POPREAL8ARRAY(cross, realtype*3/8)
   CALL CROSS_PRODUCT_3D_B(radialvec(:, i), radialvecb(:, i), &
   &                       radialvec(:, 1), radialvecb(:, 1), cross, crossb&
   &                      )
   END IF
   END DO
   CALL POPREAL8ARRAY(normal, realtype*3/8)
   centerb = 0.0_8
   DO i=npts,1,-1
   pointsb(:, i) = pointsb(:, i) + radialvecb(:, i)
   centerb = centerb - radialvecb(:, i)
   radialvecb(:, i) = 0.0_8
   END DO
   centerb = centerb/npts
   DO i=npts,1,-1
   pointsb(:, i) = pointsb(:, i) + centerb
   END DO
   END SUBROUTINE GETELEMENTPROPS_B
   