subroutine getSurfaceCoordinates(coordinates, cdof)

#include <petscversion.h>
  use gridData
  implicit none

  ! Input Arguments
  integer(kind=intType) ,  intent(in) :: cdof

  ! Output Arguments
  real(kind=realType)   ,  intent(inout) :: coordinates(cdof)

  ! Local Arguments
  integer(kind=intType) :: ierr, istart, iend, i


  call VecGetOwnershipRange(Xs, istart, iend, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  do i=1, cdof
#if PETSC_VERSION_MINOR > 13
      call VecGetValues(Xs, 1, iStart+i-1, coordinates(i), ierr)
#else
      call VecGetValues(Xs, 1, (/iStart+i-1/), coordinates(i), ierr)
#endif
    call EChk(ierr, __FILE__, __LINE__)
  end do

end subroutine getSurfaceCoordinates
