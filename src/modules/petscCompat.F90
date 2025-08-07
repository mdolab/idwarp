module petscCompat
! Include PETSc version macros
#include <petsc/finclude/petsc.h>
    ! Module for PETSc compatibility functions
    ! This module provides compatibility for various subroutines in PETSc
    ! This is because of the changes in the API between different versions of PETSc
    ! This module is used to ensure that the code works with different versions of PETSc
    use petscvec
    use constants
    implicit none

contains

    subroutine compatVecGetArray(vec, array, ierr)
        Vec :: vec
        real(kind=realType), pointer :: array(:)
        integer(kind=intType) :: ierr

#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR < 20)
        call VecGetArrayF90(vec, array, ierr)
#else
        call VecGetArray(vec, array, ierr)
#endif

    end subroutine compatVecGetArray

    subroutine compatVecRestoreArray(vec, array, ierr)
        Vec :: vec
        real(kind=realType), pointer :: array(:)
        integer(kind=intType) :: ierr

#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR < 20)
        call VecRestoreArrayF90(vec, array, ierr)
#else
        call VecRestoreArray(vec, array, ierr)
#endif

    end subroutine compatVecRestoreArray

end module petscCompat
