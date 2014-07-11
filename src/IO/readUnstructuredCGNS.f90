! ====================================================================
! File: readUnstructuredCGNS.f90
! Author: C.A.(Sandy) Mader
! Date Started: July 10, 2014
! Date Modified:

subroutine readUnstructuredCGNSFile(cgns_file, comm)

  use precision
  use communication
  implicit none
  include 'cgnslib_f.h'

  ! Input Arguments
  character*(*) :: cgns_file
  integer(kind=intType) :: comm
  
  integer(kind=intType):: ierr
  ! integer(kind=intType):: index_base, index_zone
  ! integer(kind=intType), dimension(3,1):: isize
  ! integer(kind=intType), dimension(3):: irmin, irmax

!  character*32:: zonename
  real(kind=realType), dimension(:,:,:), allocatable :: x,y,z
  ! CGNS Variabls
  integer(kind=intType) ::  cg, nzones, base, nbases

  ! ---------------------------------------
  !            Setup MPI Data first
  ! ---------------------------------------

  warp_comm_world = comm
  warp_comm_self  = mpi_comm_self
  call MPI_Comm_size(warp_comm_world, nProc, ierr)
  call MPI_Comm_rank(warp_comm_world, myid , ierr)

  ! ---------------------------------------
  !           Open CGNS File
  ! ---------------------------------------
  if (myID == 0) then
     print *, ' -> Reading CGNS File: ', cgns_file
  end if

  call cg_open_f(cgns_file, CG_MODE_READ, cg, ierr)
  if (ierr .eq. CG_ERROR) call cg_error_exit_f

  call cg_nbases_f(cg, nbases, ierr)
  if (ierr .eq. CG_ERROR) call cg_error_exit_f
  
  if (nbases .gt. 1) then
     print *, ' ** Warning: pyWarp only reads the first base in a cgns file'
  end if
  
  base = 1_intType

    !       *** base attribute:  GOTO base node
  call cg_goto_f(cg, base, ierr, 'end')
  if (ierr .eq. CG_ERROR) call cg_error_exit_f
  call cg_nzones_f(cg, base, nzones, ierr)
  if (ierr .eq. ERROR) call cg_error_exit_f
  if (myID == 0) then
     print *, '   -> Number of Zones:', nzones
  end if

! !   we know there is only one zone (real working code would check!)
!       index_zone=1
! !   get zone size (and name - although not needed here)
!       call cg_zone_read_f(index_file,index_base,index_zone,zonename,&
!            isize,ierr)
! !   lower range index
!       irmin(1)=1
!       irmin(2)=1
!       irmin(3)=1
! !   upper range index of vertices
!       irmax(1)=isize(1,1)
!       irmax(2)=isize(2,1)
!       irmax(3)=isize(3,1)
!       allocate(x(isize(1,1),isize(2,1),isize(3,1)),STAT=ierr)
!       allocate(y(isize(1,1),isize(2,1),isize(3,1)),STAT=ierr)
!       allocate(z(isize(1,1),isize(2,1),isize(3,1)),STAT=ierr)
! !   read grid coordinates
!       call cg_coord_read_f(index_file,index_base,index_zone,&
!            'CoordinateX',RealDouble,irmin,irmax,x,ierr)
!       call cg_coord_read_f(index_file,index_base,index_zone,&
!            'CoordinateY',RealDouble,irmin,irmax,y,ierr)
!       call cg_coord_read_f(index_file,index_base,index_zone,&
!            'CoordinateZ',RealDouble,irmin,irmax,z,ierr)
! !   close CGNS file
!       call cg_close_f(index_file,ierr)
!       deallocate(x,y,z,STAT=ierr)
    end subroutine readUnstructuredCGNSFile
