subroutine EChk(ierr, file, line)

  ! Check if ierr that resulted from a petsc or MPI call is in fact an
  ! error.
  use precision
#include <petscversion.h>
#if PETSC_VERSION_GE(3,8,0)
#include <petsc/finclude/petsc.h>
  use petsc
  implicit none
#else
  implicit none
#include "petsc/finclude/petsc.h"
#include "petsc/finclude/petscvec.h90"
#endif

  integer(kind=intType),intent(in) :: ierr
  character*(*),intent(in) :: file
  integer(kind=intType),intent(in) :: line
  integer(kind=intType) :: myid, jerr
  if (ierr == 0) then
     return ! No error, return immediately
  else
     call MPI_Comm_rank(petsc_comm_world, myid, jerr)
     print *,'================================================================='
     write(*,900) "PETSc or MPI Error. Error Code ",ierr,". Detected on Proc ",myid
     write(*,901) "Error at line: ",line," in file: ",file
     print *,'================================================================='

     call MPI_Abort(petsc_comm_world,ierr, jerr)
     stop ! Just in case
  end if

900 format(A,I2,A,I2)
901 format(A,I5,A,A)
end subroutine EChk
