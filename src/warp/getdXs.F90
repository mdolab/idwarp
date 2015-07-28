subroutine getdXs(indices, idof, output, ndof)

  use gridData
  implicit none
  
  ! Input Arguments
  integer(kind=intType) ,  intent(in) :: idof, ndof
  integer(kind=intType) ,  intent(in) :: indices(idof)
  
  ! Output Arguments
  real(kind=realType)   ,  intent(inout) :: output(ndof)
  
  ! Local Arguments
  integer(kind=intType) :: ierr, istart, iend, ind(idof)
  
  call VecGetOwnershipRange(dXs, istart, iend, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Pull the coordiantes out of localPatches
  output = zero
  ind = indices + istart

  call VecGetValues(dXs, ndof,  ind, output, ierr)
  call EChk(ierr, __FILE__, __LINE__)

end subroutine getdXs
