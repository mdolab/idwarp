subroutine getdXs(output, ndof)

#include <petscversion.h>
  use gridData
  implicit none

  ! Input Arguments
  integer(kind=intType) ,  intent(in) :: ndof

  ! Output Arguments
  real(kind=realType)   ,  intent(inout) :: output(ndof)

  ! Local Arguments
  integer(kind=intType) :: ierr, istart, iend, i

  call VecGetOwnershipRange(dXs, istart, iend, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  do i=1, ndof
#if PETSC_VERSION_MINOR > 13
      call VecGetValues(dXs, 1,  i+istart-1, output(i), ierr)
#else
      call VecGetValues(dXs, 1,  (/i+istart-1/), output(i), ierr)
#endif
    call EChk(ierr, __FILE__, __LINE__)
  end do

end subroutine getdXs
