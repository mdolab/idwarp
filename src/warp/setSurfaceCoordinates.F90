subroutine setSurfaceCoordinates(indices, idof, coordinates, cdof)
  
  use gridData
  implicit none
  
  ! Input Arguments
  integer(kind=intType) , intent(in) :: idof, cdof
  integer(kind=intType) , intent(in) :: indices(idof)
  real(kind=realType)   , intent(in) :: coordinates(cdof)
  
  ! Local Arguments
  integer(kind=intType) ::  ierr, istart, iend, ind(idof)
 
  call VecGetOwnershipRange(Xs, istart, iend, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  
  ind = indices + istart
  call VecSetValues(Xs, cdof, ind, coordinates, INSERT_VALUES, ierr)
  call EChk(ierr, __FILE__, __LINE__)
 
  ! While we only set local values, we STILL have to call
  ! VecAssemblyBegin/End
  call VecAssemblyBegin(Xs, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  call VecAssemblyEnd(Xs, ierr)
  call EChk(ierr, __FILE__, __LINE__)

end subroutine setSurfaceCoordinates
