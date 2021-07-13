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

  call VecGetValues(dXs, ndof, (/(i, i=istart,iend,1)/), output(1:ndof), ierr)
  call EChk(ierr, __FILE__, __LINE__)

end subroutine getdXs
