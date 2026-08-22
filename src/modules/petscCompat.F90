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

    ! Private copy of the ownership ranges handed out by VecGetOwnershipRangesCompat.
    ! PETSc >=3.23 returns a pointer that aliases the vector's internal PetscLayout and is
    ! indexed from 1, whereas older versions copied into an array the caller allocated. We
    ! always hand back a private 0:nProc copy so that callers see the same bounds on every
    ! PETSc version and can never write through to PETSc's own layout.
    integer(kind=intType), allocatable, target :: ownershipRanges(:)

contains

    subroutine VecGetArrayCompat(v, array, ierr)
        ! PETSc compatibility routine to get array from vector
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
        ! PETSc compatibility routine to restore array
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

    subroutine VecCreateMPIWithArrayCompat(comm, bs, n, nGlobal, array, v, ierr)
        ! PETSc compatibility routine to create a MPI vector with an array.
        ! In some cases array is not supplied so it supports optional array argument.
        implicit none
        integer(kind=intType), intent(in) :: comm, bs
        integer(kind=intType), intent(in) :: n, nGlobal
        real(kind=realType), dimension(:), intent(in), optional :: array
        Vec :: v
        integer(kind=intType) :: ierr

#if PETSC_VERSION_GE(3,22,0)
        ! PETSc >= 3.22, need to use PETSC_NULL_SCALAR_ARRAY if array is not present
        if (present(array)) then
            call VecCreateMPIWithArray(comm, bs, n, nGlobal, array, v, ierr)
        else
            call VecCreateMPIWithArray(comm, bs, n, nGlobal, PETSC_NULL_SCALAR_ARRAY, v, ierr)
        end if
#else
        ! Older PETSc, use PETSC_NULL_SCALAR if array is not present
        if (present(array)) then
            call VecCreateMPIWithArray(comm, bs, n, nGlobal, array, v, ierr)
        else
            call VecCreateMPIWithArray(comm, bs, n, nGlobal, PETSC_NULL_SCALAR, v, ierr)
        end if
#endif

    end subroutine VecCreateMPIWithArrayCompat

    subroutine VecGetOwnershipRangesCompat(v, ptr, ierr)
        ! PETSc compatibility routine to get ownership ranges.
        ! On every PETSc version this returns a pointer to a private copy indexed 0:nProc,
        ! so that ptr(iProc) is the first index owned by rank iProc and ptr(nProc) is the
        ! global size. The caller may modify the copy freely and must release it with
        ! VecRestoreOwnershipRangesCompat.
        Vec :: v
        integer(kind=intType), pointer :: ptr(:)
        integer(kind=intType), intent(out) :: ierr

#if PETSC_VERSION_GE(3,23,0)
        integer(kind=intType), pointer :: petscRanges(:)
#endif

        ! Discard any copy left over from a previous call
        if (allocated(ownershipRanges)) then
            deallocate (ownershipRanges)
        end if
        ! Use zero-based indexing to match to rank
        allocate (ownershipRanges(0:nProc))

#if PETSC_VERSION_GE(3,23,0)
        ! PETSc >= 3.23 hands back a pointer aliasing its own PetscLayout, indexed from 1.
        ! Copy it out and give it straight back, rebasing to 0:nProc as we go. Writing
        ! through the PETSc pointer would silently corrupt the vector's layout.
        call VecGetOwnershipRanges(v, petscRanges, ierr)
        if (ierr == 0) then
            ownershipRanges = petscRanges
            call VecRestoreOwnershipRanges(v, petscRanges, ierr)
        end if
#else
        ! Older PETSc copies into an array we own, so hand it ours directly
        call VecGetOwnershipRanges(v, ownershipRanges, ierr)
#endif

        ptr => ownershipRanges

    end subroutine VecGetOwnershipRangesCompat

    subroutine VecRestoreOwnershipRangesCompat(v, ptr, ierr)
        ! PETSc compatibility routine to release the ownership ranges copy.
        ! The PETSc-owned pointer was already returned inside VecGetOwnershipRangesCompat,
        ! so on all versions this just frees our private copy.
        Vec :: v
        integer(kind=intType), pointer :: ptr(:)
        integer(kind=intType), intent(out) :: ierr

        if (allocated(ownershipRanges)) then
            deallocate (ownershipRanges)
        end if
        nullify (ptr)
        ierr = 0

    end subroutine VecRestoreOwnershipRangesCompat

end module petscCompat
