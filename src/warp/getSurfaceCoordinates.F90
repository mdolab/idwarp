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

  do i=1, cdof
     call VecGetValues(Xs, 1, iStart+i-1, coordinates(i), ierr)
     call EChk(ierr, __FILE__, __LINE__)
  end do

end subroutine getSurfaceCoordinates
