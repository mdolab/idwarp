subroutine setdXs(indices, idof, in_vec, ndof)

  use precision 
  use gridData
  implicit none
  
  ! Input Arguments
  integer(kind=intType) , intent(in) :: ndof, idof
  integer(kind=intType) , intent(in) :: indices(idof)
  real(kind=realType)   , intent(in) :: in_vec(ndof)

  ! Local Arguments
  integer(kind=intType) :: ierr, istart, iend, ind(idof)

  call VecGetOwnershipRange(dxs, istart, iend, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Zero dXs
  call VecZeroEntries(dXs, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  
  ind = istart + indices
  ! Set provided values at providedindices
  call VecSetValues(dxs, ndof, ind, in_vec, INSERT_VALUES, ierr)
  call EChk(ierr, __FILE__, __LINE__)
 
  ! While we only set local values, we STILL have to call
  ! VecAssemblyBegin/End
  call VecAssemblyBegin(dxs, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  call VecAssemblyEnd(dxs, ierr)
  call EChk(ierr, __FILE__, __LINE__)

end subroutine setdXs

