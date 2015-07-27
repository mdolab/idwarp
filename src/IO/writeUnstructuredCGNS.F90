subroutine writeUnstructuredCGNS(cgns_file)

  use precision
  use communication
  use unstructuredcgnsgrid
  use gridData
  implicit none
  include 'cgnslib_f.h'

  ! Input Arguments
  character*(*),intent(in) :: cgns_file
    
  ! CGNS Variabls
  integer(kind=intType) :: cg, base
  character*32 :: zoneName, bocoName, famName, connectName, donorName, baseName
  integer(kind=intType) :: iZone, zone, sizes(3), coordID, ierr

  real(kind=realType), dimension(:), allocatable :: coorX, coorY, coorZ
  real(kind=realType), pointer, dimension(:) :: xx
  integer(kind=intType) :: i, j

  ! This isn't technically scalable...we will dump the entire grid
  ! onto the root proc and write there. 
  call VecScatterCreateToZero(Xv, XvToLocal, Xvlocal, ierr)

  call VecScatterBegin(XvToLocal, Xv, XvLocal, INSERT_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecScatterEnd(XvToLocal, Xv, XvLocal, INSERT_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Only do writing on root proc:
  rootProc: if (myid == 0) then 

     call VecGetArrayF90(XvLocal, xx, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Open and get the number of zones:
     call cg_open_f(trim(cgns_file), CG_MODE_MODIFY, cg, ierr)
     if (ierr .eq. CG_ERROR) call cg_error_exit_f
     base = 1

     ! loop over the zones, create and write
     do iZone=1, nZones

        ! Get the number of points and elements
        sizes(1) = zones(iZone)%nVertices
        sizes(2) = zones(iZone)%nElements
        sizes(3) = 0

        ! Allocate arrays for X,Y and Z.Note that we CANNOT just use
        ! matrix notation and stride directly from xx since that will
        ! cause a stack copy to be creater which will make ifort barf,
        ! since they put temps are the stack by default *shakes fist*
        allocate(coorX(sizes(1)), coorY(sizes(1)), coorZ(sizes(1)))

        ! Extract the coordinates
        do j=1,sizes(1)
           coorX(j) = xx(3*j-2)
           coorY(j) = xx(3*j-1)
           coorZ(j) = xx(3*j  )
        end do

        ! Write the grid coordinates
        call cg_coord_write_f(cg, base, iZone, realDouble, 'CoordinateX', coorX, coordID, ierr)
        if (ierr .eq. CG_ERROR) call cg_error_exit_f

        call cg_coord_write_f(cg, base, iZone, realDouble, 'CoordinateY', coorY, coordID, ierr)
        if (ierr .eq. CG_ERROR) call cg_error_exit_f

        call cg_coord_write_f(cg, base, iZone, realDouble, 'CoordinateZ', coorZ, coordID, ierr)
        if (ierr .eq. CG_ERROR) call cg_error_exit_f

        deallocate(coorX, coorY, coorZ)
     end do
     call cg_close_f(cg, ierr)
     if (ierr .eq. CG_ERROR) call cg_error_exit_f

     call VecRestoreArrayF90(XvLocal, xx, ierr)
     call EChk(ierr,__FILE__,__LINE__)

  end if rootProc
  call vecDestroy(XvLocal, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecScatterDestroy(XvToLocal, ierr)
  call EChk(ierr,__FILE__,__LINE__)

end subroutine writeUnstructuredCGNS

