subroutine warpDerivFwd(indices, idof, Xsdot, cDof, Xvdot, meshDOF)

  use gridData
  use gridInput
  use communication
  use kd_tree
  implicit none
  
  ! Input Arguments
  integer(kind=intType) , intent(in) :: idof, cdof, meshDOF
  integer(kind=intType) , intent(in) :: indices(idof)
  real(kind=realType)   , intent(in) :: xsdot(cdof)
  
  ! Output Arguments
  real(kind=realType), intent(inout), dimension(meshDOF) :: XvDot
#ifndef USE_COMPLEX
  ! Working Data
  integer(kind=intType) :: i, j, kk, istart, iend, ierr, isize, nVol, nLoop
  integer(kind=intType) :: ind(idof)
  real(kind=realType), dimension(3) :: r, rr, numd
  real(kind=realType) :: fact(3, 2), oden

  ! Dump our Xsdot into the dXs array
  call VecZeroEntries(dXs, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecGetOwnershipRange(dXs, istart, iend, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  
  ind = indices + istart
  call VecSetValues(dXs, cdof, ind, XsDot, INSERT_VALUES, ierr)
  call EChk(ierr, __FILE__, __LINE__)
 
  ! While we only set local values, we STILL have to call
  ! VecAssemblyBegin/End
  call VecAssemblyBegin(dXs, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  call VecAssemblyEnd(dXs, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Scatter Xs into our local vector  
  call VecScatterBegin(XsToXsLocal, dXs, dXsLocal, INSERT_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecScatterEnd(XsToXsLocal, dXs, dXsLocal, INSERT_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Extract a pointer from XsLocal and XsLoacld to use in the main routine
  call VecGetArrayF90(XsLocal, XsPtr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecGetArrayF90(dXsLocal, XsPtrd, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Extract the pointers for the volume nodes that we wish to operate on
  call VecGetArrayF90(Xv0, Xv0Ptr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call allocDerivValues(mytrees(1)%tp)
  call zeroDeriv(mytrees(1)%tp%root)

  nVol = size(XvPtr)/3
  call computeNodalProperties_d(mytrees(1)%tp, .False.)
  mytrees(1)%tp%Ldef = mytrees(1)%tp%Ldef0 * LdefFact
  mytrees(1)%tp%alphaToBexp = alpha**bExp
  mytrees(1)%tp%errTol = errTol
  call setData_d(mytrees(1)%tp, mytrees(1)%tp%root)

  call getLoopFact(nLoop, fact)

  numerator = zero
  do kk=1, nLoop
     volLoop: do j=1, nVol
        r(1) = Xv0Ptr(3*j-2)*fact(1, kk)
        r(2) = Xv0Ptr(3*j-1)*fact(2, kk)
        r(3) = Xv0Ptr(3*j  )*fact(3, kk)

        numd = zero
        if (evalMode == EVAL_EXACT) then 
           call evalNodeExact_d(mytrees(1)%tp, mytrees(1)%tp%root, r, numd)
        else           
           call evalNode_d(mytrees(1)%tp, mytrees(1)%tp%root, r, numd, &
                denomenator0(j))
        end if
        ! Sum the *derivatives* into numerator
        numerator(1, j) = numerator(1, j) + numd(1)*fact(1, kk)
        numerator(2, j) = numerator(2, j) + numd(2)*fact(2, kk)
        numerator(3, j) = numerator(3, j) + numd(3)*fact(3, kk)

     end do volLoop
  end do
  
  ! Done with derivative values 
  call deallocDerivValues(mytrees(1)%tp)
  
  updateLoop: do j=1, nVol
     oden = one / denomenator(j)
     Xvdot(3*j-2) = numerator(1, j) * oden
     Xvdot(3*j-1) = numerator(2, j) * oden
     Xvdot(3*j  ) = numerator(3, j) * oden
  end do updateLoop

  ! Restore all the arrays
  call VecRestoreArrayF90(XsLocal, XsPtr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecRestoreArrayF90(dXsLocal, dXsPtr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecRestoreArrayF90(Xv0, Xv0Ptr, ierr)
  call EChk(ierr,__FILE__,__LINE__)
#endif
end subroutine warpDerivFwd
