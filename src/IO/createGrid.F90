subroutine createVolumeGrid(volNodes, nVolLocal)

  use gridData
  use communication
  implicit none

  ! Subroutine variables
  real(kind=realType), dimension(3, nVolLocal), intent(in) :: volNodes
  integer(kind=intType), intent(in) :: nVolLocal

  ! Working variables
  integer(kind=intType) :: i, nVol, ierr
  integer(kind=intType), dimension(:), allocatable :: volNodesProc

  ! Gather the displacements
  allocate(volNodesProc(0:nProc))
  volNodesProc(:) = 0
  call MPI_allgather(nVolLocal, 1, MPI_INTEGER, volNodesProc(1:), 1, MPI_INTEGER, &
       warp_comm_world, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Finish the displacment calc:
  do i=2, nProc
     volNodesProc(i) = volNodesProc(i) + volNodesProc(i-1) 
  end do
  nVol = volNodesProc(nProc)

  if (myid == 0) then
     write(*,"(a)", advance="no") '#------------------------------#'
     print "(1x)"  
     write(*,"(a)", advance="no") " Total Volume Nodes : "
     write(*,"(I9,1x)",advance="no") nVol
     print "(1x)"  
     write(*,"(a)", advance="no") '#------------------------------#'
     print "(1x)"   
  end if

  ! Now create the master Xv array:
  call VecCreate(warp_comm_world, Xv, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Set to be be blocked
  call VecSetBlockSize(Xv, 3, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Type and size
  call VecSetType(Xv, "mpi", ierr)
  call EChk(ierr,__FILE__,__LINE__)
  
  call VecSetSizes(Xv, nVolLocal*3, PETSC_DECIDE, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Now each processor adds it's own nodes (It is possible only one proc does it)
  do i=1, nVolLocal 
     call VecSetValuesBlocked(Xv, 1, volNodesProc(myid) + i - 1, volNodes(:, i), INSERT_VALUES, ierr)
     call EChk(ierr,__FILE__,__LINE__)
  end do

  call VecGetLocalSize(Xv, warpMeshDOF, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecAssemblyBegin(Xv, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecAssemblyEnd(Xv, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecDuplicate(Xv, Xv0, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecCopy(Xv, Xv0, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  deallocate(volNodesProc)
  allocate(numerator(3, nVolLocal))
  allocate(denomenator(nVolLocal), denomenator0(nVolLocal))

end subroutine createVolumeGrid


