subroutine releaseMemory
  ! Clean up all fortran-allocated memory. 
  use kd_tree
  use gridData
  use constants
  use CGNSGrid
  implicit none

  ! Working parameters
  integer(kind=intType) :: i, j, ierr

  ! Created in 'createCommonGrid'
  if (commonGridVecSet == 1) then 
     call VecDestroy(commonGridVec, ierr)
     call EChk(ierr,__FILE__,__LINE__)

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

  ! Check if we zones from a CGNS mesh
  if (allocated(zones)) then 
     do i=1,size(zones)
        
        ! We have bocos
        if (allocated(zones(i)%bocos)) then 
           
           ! Loop over all the bocos
           do j=1,size(zones(i)%bocos)

              ! Delete allocatable data in Boco
              if (allocated(zones(i)%bocos(j)%BCElements)) then
                 deallocate(zones(i)%bocos(j)%BCElements)
              end if

              if (associated(zones(i)%bocos(j)%elemPtr)) then 
                 deallocate(zones(i)%bocos(j)%elemPtr)
              end if

              if (associated(zones(i)%bocos(j)%elemNodes)) then 
                 deallocate(zones(i)%bocos(j)%elemNodes)
              end if

              if (associated(zones(i)%bocos(j)%elemConn)) then 
                 deallocate(zones(i)%bocos(j)%elemConn)
              end if
           end do

           ! Delete boco array itself
           deallocate(zones(i)%bocos)
        end if
        
        ! Sections also take more work
        if (associated(zones(i)%sections)) then 

           ! Loop over all sections
           do j=1, size(zones(i)%sections)
              if (associated(zones(i)%sections(j)%elemPtr)) then 
                 deallocate(zones(i)%sections(j)%elemPtr)
              end if

              if (associated(zones(i)%sections(j)%elemConn)) then 
                 deallocate(zones(i)%sections(j)%elemConn)
              end if
           end do
           
           ! Delete section array itself
           deallocate(zones(i)%sections)
        end if
     end do

     ! Delete zone array itself
     deallocate(zones)

     ! Delete the surface-stuff in the cgns grid if necessary
     if (allocated(surfacePoints)) then 
        deallocate(surfacePoints)
     end if

     if (allocated(surfaceConn)) then 
        deallocate(surfaceConn)
     end if

     if (allocated(surfaceSizes)) then 
        deallocate(surfaceSizes)
     end if

     if (allocated(surfacePtr)) then 
        deallocate(surfacePtr)
     end if

     if (allocated(surfacePatchPtr)) then 
        deallocate(surfacePatchPtr)
     end if

     if (allocated(surfaceNames)) then 
        deallocate(surfaceNames)
     end if

     if (allocated(surfaceIsWall)) then 
        deallocate(surfaceIsWall)
     end if

     if (allocated(surfaceIsSymm)) then 
        deallocate(surfaceIsSymm)
     end if

  end if ! Using CGNS grid

  if (allocated(symmPts)) then 
     deallocate(symmPts)
  end if
  
  if (allocated(symmNormals)) then 
     deallocate(symmNormals)
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
