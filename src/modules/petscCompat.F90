module petscCompat
! Include PETSc version macros
#include <petsc/finclude/petsc.h>
    ! Module for PETSc compatibility functions
    ! This module provides compatibility for various subroutines in PETSc
    ! This is because of the changes in the API between different versions of PETSc
    ! This module is used to ensure that the code works with different versions of PETSc

    use petscvec
    use constants
    use communication
    implicit none
    save

    ! Internal static array needed for PETSc <3.23 VecGetOwnershipRanges
    integer(kind=intType), allocatable, target :: array(:)
    logical :: arrayAllocated = .false.

contains

    subroutine VecGetArrayCompat(v, array, ierr)
        ! PETSc compatability routine to get array from vector
        implicit none
        Vec :: v
        real(kind=realType), pointer :: array(:)
        integer(kind=intType) :: ierr

#if PETSC_VERSION_LT(3,20,0)
        call VecGetArrayF90(v, array, ierr)
#else
        call VecGetArray(v, array, ierr)
#endif

    end subroutine VecGetArrayCompat

    subroutine VecRestoreArrayCompat(v, array, ierr)
        ! PETSc compatability routine to restore array
        implicit none
        Vec :: v
        real(kind=realType), pointer :: array(:)
        integer(kind=intType) :: ierr

#if PETSC_VERSION_LT(3,20,0)
        call VecRestoreArrayF90(v, array, ierr)
#else
        call VecRestoreArray(v, array, ierr)
#endif

    end subroutine VecRestoreArrayCompat

    subroutine VecCreateMPIWithArrayCompat(comm, bs, nlocal, nglobal, array, v, ierr)
        ! PETSc compatability routine to create a MPI vector with an array.
        ! In some cases array is not supplied so it supports optional array argument.
        implicit none
        integer(kind=intType), intent(in) :: comm, bs
        integer(kind=intType), intent(in) :: nlocal, nglobal
        real(kind=realType), dimension(:), intent(in), optional :: array
        Vec :: v
        integer(kind=intType) :: ierr

#if PETSC_VERSION_GE(3,22,0)
        ! PETSc >= 3.22, need to use PETSC_NULL_SCALAR_ARRAY
        if (present(array)) then
            call VecCreateMPIWithArray(comm, bs, nlocal, nglobal, array, v, ierr)
        else
            call VecCreateMPIWithArray(comm, bs, nlocal, nglobal, PETSC_NULL_SCALAR_ARRAY, v, ierr)
        end if
#else
        ! Older PETSc, use PETSC_NULL_SCALAR for older versions
        if (present(array)) then
            call VecCreateMPIWithArray(comm, bs, nlocal, nglobal, array, v, ierr)
        else
            call VecCreateMPIWithArray(comm, bs, nlocal, nglobal, PETSC_NULL_SCALAR, v, ierr)
        end if
#endif

    end subroutine VecCreateMPIWithArrayCompat

    subroutine VecGetOwnershipRangesCompat(v, ptr, ierr)
        ! PETSc compatability routine to get ownership ranges
        Vec :: v
        integer(kind=intType), pointer :: ptr(:)
        integer(kind=intType), intent(out) :: ierr

#if PETSC_VERSION_GE(3,23,0)
        ! PETSc >= 3.23, directly get pointer
        call VecGetOwnershipRanges(v, ptr, ierr)
#else
        ! Older PETSc, allocate array that gets passed to PETSc and set pointer
        integer(kind=intType), allocatable, target, save :: array(:)

        ! Check if we already have an allocated array
        if (allocated(array)) then
            deallocate (array)
        end if
        ! Use zero-based indexing to match to rank
        allocate (array(0:nProc))
        array = zero
        arrayAllocated = .true.

        ! Get the ownership ranges into the static array and set the pointer
        call VecGetOwnershipRanges(v, array, ierr)
        ptr => array
#endif

    end subroutine VecGetOwnershipRangesCompat

    subroutine VecRestoreOwnershipRangesCompat(v, ptr, ierr)
        ! PETSc compatability routine to restore ownership ranges pointer
        Vec :: v
        integer(kind=intType), pointer :: ptr(:)
        integer(kind=intType), intent(out) :: ierr

#if PETSC_VERSION_GE(3,23,0)
        ! PETSc >= 3.23, just restore ranges using the pointer
        call VecRestoreOwnershipRanges(v, ptr, ierr)
#else
        ! Older PETSc, deallocate the internal static array and nullify pointer
        if (arrayAllocated) then
            if (allocated(array)) then
                deallocate (array)
            end if
            nullify (ptr)
            arrayAllocated = .false.
        end if
        ierr = 0
#endif

    end subroutine VecRestoreOwnershipRangesCompat

end module petscCompat
