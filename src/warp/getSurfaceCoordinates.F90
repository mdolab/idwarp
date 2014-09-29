subroutine getSurfaceCoordinates(indices, idof, coordinates, cdof)

  use gridData
  implicit none
  
  ! Input Arguments
  integer(kind=intType) ,  intent(in) :: idof, cdof
  integer(kind=intType) ,  intent(in) :: indices(idof)
  
  ! Output Arguments
  real(kind=realType)   ,  intent(inout) :: coordinates(cdof)
  
  ! Local Arguments
  integer(kind=intType) :: ierr, istart, iend, ind(idof)
  
  call VecGetOwnershipRange(Xs, istart, iend, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Pull the coordiantes out of localPatches
  coordinates = zero

  ind = indices + istart
  call VecGetValues(Xs, idof,  ind, coordinates, ierr)
  call EChk(ierr, __FILE__, __LINE__)

end subroutine getSurfaceCoordinates
