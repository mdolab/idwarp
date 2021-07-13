subroutine getSurfaceCoordinates(coordinates, cdof)

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

  call VecGetValues(Xs, cdof, (/(i, i=istart,iend,1)/), coordinates(1:cdof), ierr)
  call EChk(ierr, __FILE__, __LINE__)

end subroutine getSurfaceCoordinates
