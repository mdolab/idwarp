subroutine warpDerivFwd(Xsdot, cDof, Xvdot, meshDOF)

    use gridData
    use gridInput
    use kd_tree

    implicit none

    ! Input Arguments
    integer(kind=intType), intent(in) :: cdof, meshDOF
    real(kind=realType), intent(in) :: xsdot(cdof)

    ! Output Arguments
    real(kind=realType), intent(inout), dimension(meshDOF) :: XvDot
#ifndef USE_COMPLEX
    ! Working Data
    integer(kind=intType) :: i, j, kk, istart, iend, ierr, nVol
    real(kind=realType), dimension(3) :: r, numd
    real(kind=realType) :: oden

    ! Dump our Xsdot into the dXs array
    call VecZeroEntries(dXs, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call VecGetOwnershipRange(dXs, istart, iend, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    do i = 1, cdof
        call VecSetValues(dXs, 1, [iStart + i - 1], [XsDot(i)], INSERT_VALUES, ierr)
        call EChk(ierr, __FILE__, __LINE__)
    end do

    ! While we only set local values, we STILL have to call
    ! VecAssemblyBegin/End
    call VecAssemblyBegin(dXs, ierr)
    call EChk(ierr, __FILE__, __LINE__)
    call VecAssemblyEnd(dXs, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Scatter Xs into our local vector
    call VecScatterBegin(XsToXsLocal, dXs, dXsLocal, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call VecScatterEnd(XsToXsLocal, dXs, dXsLocal, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Extract a pointer from XsLocal and XsLoacld to use in the main routine
    call VecGetArrayF90(XsLocal, XsPtr, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call VecGetArrayF90(dXsLocal, XsPtrd, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Extract the pointers for the volume nodes that we wish to operate on
    call VecGetArrayF90(Xv0, Xv0Ptr, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call allocDerivValues(mytrees(1)%tp)
    call zeroDeriv(mytrees(1)%tp%root)

    nVol = size(XvPtr) / 3
    call computeNodalProperties_d(mytrees(1)%tp, .False.)
    mytrees(1)%tp%Ldef = mytrees(1)%tp%Ldef0 * LdefFact
    mytrees(1)%tp%alphaToBexp = alpha**bExp
    mytrees(1)%tp%errTol = errTol
    call setData_d(mytrees(1)%tp, mytrees(1)%tp%root)

    numerator = zero
    do kk = 1, nLoop
        volLoop: do j = 1, nVol
            call getMirrorPt(Xv0Ptr(3 * j - 2:3 * j), r, kk)

            numd = zero
            if (evalMode == EVAL_EXACT) then
                call evalNodeExact_d(mytrees(1)%tp, mytrees(1)%tp%root, r, numd)
            else
                call evalNode_d(mytrees(1)%tp, mytrees(1)%tp%root, r, numd, &
                                denomenator0(j))
            end if

            call getMirrorNumerator_d(numd, kk)

            ! Sum the *derivatives* into numerator
            numerator(1, j) = numerator(1, j) + numd(1)
            numerator(2, j) = numerator(2, j) + numd(2)
            numerator(3, j) = numerator(3, j) + numd(3)

        end do volLoop
    end do

    ! Done with derivative values
    call deallocDerivValues(mytrees(1)%tp)

    updateLoop: do j = 1, nVol
        oden = one / denomenator(j)
        Xvdot(3 * j - 2) = numerator(1, j) * oden
        Xvdot(3 * j - 1) = numerator(2, j) * oden
        Xvdot(3 * j) = numerator(3, j) * oden
    end do updateLoop

    ! Restore all the arrays
    call VecRestoreArrayF90(XsLocal, XsPtr, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call VecRestoreArrayF90(dXsLocal, dXsPtr, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call VecRestoreArrayF90(Xv0, Xv0Ptr, ierr)
    call EChk(ierr, __FILE__, __LINE__)
#endif
end subroutine warpDerivFwd
