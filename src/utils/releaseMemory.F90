subroutine releaseMemory
  ! Clean up all fortran-allocated memory. 
  use kd_tree
  use gridData
  use constants
  use structuredCGNSGrid
  implicit none

  ! Working parameters
  integer(kind=intType) :: i, ierr

  ! Created in 'createCommonGrid'
  if (commonGridVecSet == 1) then 
     call VecDestroy(commonGridVec, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     deallocate(wallIndices)
     commonGridVecSet = 0
  end if
  
  if (initializationSet == 1) then 
     ! The following were created in initialzieWarping:
     call VecDestroy(Xs, ierr)
     call EChk(ierr,__FILE__,__LINE__)
     
     call VecDestroy(dXs, ierr)
     call EChk(ierr,__FILE__,__LINE__)
     
     call VecDestroy(XsLocal, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call VecDestroy(dXsLocal, ierr)
     call EChk(ierr,__FILE__,__LINE__)
     
     call VecScatterDestroy(XsToXsLocal, ierr)
     call EChk(ierr,__FILE__,__LINE__)
     
     call VecDestroy(Xv, ierr)
     call EChk(ierr,__FILE__,__LINE__)
     
     call VecDestroy(Xv0, ierr)
     call EChk(ierr,__FILE__,__LINE__)
     
     call VecDestroy(dXv, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call VecScatterDestroy(common_to_warp, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call VecScatterDestroy(common_to_dXs, ierr)
     call EChk(ierr,__FILE__,__LINE__)
    

     ! Deallocate the num/den data
     deallocate(numerator, denomenator, denomenator0)
     
     initializationSet = 0
  end if

  ! Check if we have a structured mesh
  if (allocated(blocks)) then 
     do i=1,size(blocks)
        if (allocated(blocks(i)%bocos)) then 
           deallocate(blocks(i)%bocos)
        end if
        if (allocated(blocks(i)%B2Bs)) then 
           deallocate(blocks(i)%B2Bs)
        end if
     end do
     deallocate(blocks)
  end if
  
  if (gridIndicesSet == 1) then
     call VecDestroy(solverGridVec, ierr)
     call EChk(ierr,__FILE__,__LINE__)
     
     call VecScatterDestroy(common_to_solver, ierr)
     call EChk(ierr,__FILE__,__LINE__)
     gridIndicesSet = 0
  end if

  ! Finally delete the trees
  if (allocated(mytrees)) then
     do i=1,size(mytrees)
        call destroy_tree(mytrees(i)%tp)
     end do
     deallocate(mytrees)
  end if
end subroutine releaseMemory
