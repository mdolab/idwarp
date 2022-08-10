subroutine verifyWarpDeriv(dXv_f, ndof_warp, dof_start, dof_end, h)

#include <petscversion.h>
    use gridData
    use gridInput
    use communication

    implicit none

    ! Input
    integer(kind=intType) :: ndof_warp
    real(kind=realType) :: dXv_f(ndof_warp), h
    integer(kind=inttype) :: dof_start, dof_end

    ! Working
    integer(kind=intType) :: istart, iend, dof, dofSurfMax, ierr, nDof_to_check
    real(kind=realType) :: FDvalue, ADValue(1), value, err, orig_value(1), val
    real(kind=realType), dimension(:), allocatable :: deriv, xplus, xminus

    allocate (deriv(warpMeshDOF), xplus(warpMeshDOF), xminus(warpMeshDOF))

    call VecGetSize(dXs, dofSurfMax, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    dof_start = max(0, min(dof_start, dofSurfMax - 1))
    dof_end = min(dofSurfMax - 1, max(dof_end, 0))

    NDof_to_check = dof_end - dof_start

    if (myid == 0) then
        print *, 'Welcome to verifyWarpDeriv'
        print *, 'Checking ', NDof_to_check, ' degrees of freedom'
        print *, 'Doing centered differnce with h:', h
    end if
    ! ---------------------------------------------
    !               AD Check
    ! ---------------------------------------------

    if (myid == 0) then
        print *, 'Running AD Version'
    end if

    call warpMesh()
    call WarpDeriv(dXv_f, ndof_warp)

    ! dXs now contains actual AD derivative
    call VecGetOwnershipRange(dXs, istart, iend, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Loop over desired DOFs
    do dof = dof_start, dof_end

        ! add h to dof
        if (dof >= istart .and. dof < iend) then
#if PETSC_VERSION_GE(3,14,0)
            call VecGetValues(Xs, 1, dof, orig_value, ierr)
#else
            call VecGetValues(Xs, 1, (/dof/), orig_value, ierr)
#endif
            call EChk(ierr, __FILE__, __LINE__)

            val = orig_value(1) + h
            call VecSetValue(Xs, (/dof/), val, INSERT_VALUES, ierr)
            call EChk(ierr, __FILE__, __LINE__)
        end if

        call VecAssemblyBegin(Xs, ierr)
        call EChk(ierr, __FILE__, __LINE__)
        call VecAssemblyEnd(Xs, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call warpMesh()

        ! Copy what is is Xv into Xplus
        call VecGetArrayF90(Xv, XvPtr, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        xplus = XvPtr

        call VecRestoreArrayF90(Xv, XvPtr, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! Subtract 2h from dof to get x(dof)-h
        if (dof >= istart .and. dof < iend) then
            val = orig_value(1) - h
            call VecSetValue(Xs, (/dof/), val, INSERT_VALUES, ierr)
            call EChk(ierr, __FILE__, __LINE__)
        end if

        call VecAssemblyBegin(Xs, ierr)
        call EChk(ierr, __FILE__, __LINE__)
        call VecAssemblyEnd(Xs, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! Warp Mesh (again)
        call warpMesh()

        ! Copy what is is Xv into Xminus
        call VecGetArrayF90(Xv, XvPtr, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        xminus = XvPtr

        call VecRestoreArrayF90(Xv, XvPtr, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! reset the original value
        if (dof >= istart .and. dof < iend) then
            call VecSetValue(Xs, (/dof/), orig_value(1), INSERT_VALUES, ierr)
            call EChk(ierr, __FILE__, __LINE__)
        end if

        call VecAssemblyBegin(Xs, ierr)
        call EChk(ierr, __FILE__, __LINE__)
        call VecAssemblyEnd(Xs, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! Do actual central differecing fd calc
        deriv = (one / (two * h)) * (xplus - xminus) ! Overwrite into deriv

        ! Now sum up which is the same as taking the dot-product with
        ! a set of ones

        value = dot_product(dxv_f, deriv)

        call MPI_allReduce(value, FDValue, 1, MPI_REAL8, MPI_SUM, &
                           warp_comm_world, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        if (dof >= istart .and. dof < iend) then
#if PETSC_VERSION_GE(3,14,0)
            call VecGetValues(dXs, 1, dof, ADvalue, ierr)
#else
            call VecGetValues(dXs, 1, (/dof/), ADvalue, ierr)
#endif
            call EChk(ierr, __FILE__, __LINE__)
            if (abs(half * (FDValue + ADValue(1))) < 1e-16) then
                err = 1e-16
            else
                err = (FDValue - ADValue(1)) / (half * (FDValue + ADValue(1))) * 100_realType
            end if

            write (*, 900) 'DOF:', dof, ' OrigVal: ', orig_value(1), ' AD:', ADValue, ' FD:', &
                FDValue, ' Err(%):', err
        end if

900     format(A, I5, A, G19.12, A, G19.12, A, G19.12, A, G17.10)
    end do

    ! Make sure mesh is up to date
    call warpMesh()

    deallocate (deriv, xplus, xminus)
end subroutine verifyWarpDeriv
