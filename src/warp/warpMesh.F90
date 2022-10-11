subroutine warpMesh()

    ! This is the master mesh warping gateway function. It does not do
    ! any actual computations, rather just extracts the required
    ! pointers from petsc, determines which of the two gateway functions
    ! are called depending on if the fast or exact code are called.

    use gridData
    use gridInput
    use kd_tree

    implicit none

    ! Working parameters
    real(kind=realType) :: den, oden
    real(kind=realType), dimension(3) :: r, num
    integer(kind=intType) :: nVol, ierr, j, kk

    ! Scatter Xs into our local vector
    call VecScatterBegin(XsToXsLocal, Xs, XsLocal, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call VecScatterEnd(XsToXsLocal, Xs, XsLocal, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Extract a pointer from Xs to use in the main routine
    call VecGetArrayF90(XsLocal, XsPtr, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Extract the pointers for the volume nodes that we wish to operate on
    call VecGetArrayF90(Xv0, Xv0Ptr, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call VecGetArrayF90(Xv, XvPtr, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    nVol = size(XvPtr) / 3
    call computeNodalProperties(mytrees(1)%tp, .False.)
    mytrees(1)%tp%Ldef = mytrees(1)%tp%Ldef0 * LdefFact
    mytrees(1)%tp%alphaToBexp = alpha**bExp
    mytrees(1)%tp%errTol = errTol
    call setData(mytrees(1)%tp, mytrees(1)%tp%root)

    denominator = zero
    numerator = zero
    do kk = 1, nLoop
        volLoop: do j = 1, nVol
            call getMirrorPt(Xv0Ptr(3 * j - 2:3 * j), r, kk)
            num = zero
            den = zero
            if (evalMode == EVAL_EXACT) then
                call evalNodeExact(mytrees(1)%tp, mytrees(1)%tp%root, r, num, den)
            else
                call evalNode(mytrees(1)%tp, mytrees(1)%tp%root, r, num, den, denominator0(j))
            end if

            call getMirrorNumerator(num, kk)
            numerator(1, j) = numerator(1, j) + num(1)
            numerator(2, j) = numerator(2, j) + num(2)
            numerator(3, j) = numerator(3, j) + num(3)
            denominator(j) = denominator(j) + den
        end do volLoop
    end do

    updateLoop: do j = 1, nVol
        oden = one / denominator(j)
        XvPtr(3 * j - 2) = Xv0Ptr(3 * j - 2) + numerator(1, j) * oden
        XvPtr(3 * j - 1) = Xv0Ptr(3 * j - 1) + numerator(2, j) * oden
        XvPtr(3 * j) = Xv0Ptr(3 * j) + numerator(3, j) * oden
    end do updateLoop

    ! Restore all the arrays
    call VecRestoreArrayF90(XsLocal, XsPtr, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call VecRestoreArrayF90(Xv0, Xv0Ptr, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call VecRestoreArrayF90(Xv, XvPtr, ierr)
    call EChk(ierr, __FILE__, __LINE__)

end subroutine warpMesh

