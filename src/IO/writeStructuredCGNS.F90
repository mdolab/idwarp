subroutine writeStructuredCGNS(cgns_file)
  use gridData
  use communication
  use CGNSGrid

  implicit none
  include 'cgnslib_f.h'

  ! Input Arguments
  character*(*),intent(in) :: cgns_file

  ! Working 
  integer(kind=intType) :: i,ii, jj
  integer(kind=intType):: ierr, base, cg, curSize, nCon, bcOut, boco, ib2b, zoneid
  integer(kind=intType) :: sizes(9), coordID
  real(kind=realType), dimension(:), allocatable :: coorX, coorY, coorZ
  real(kind=realType), pointer, dimension(:) :: xx

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
     call cg_open_f(trim(cgns_file), CG_MODE_WRITE, cg, ierr)
     if (ierr .eq. CG_ERROR) call cg_error_exit_f

     call cg_base_write_f(cg, "BASE#1", 3, 3, base, ierr)
     if (ierr .eq. CG_ERROR) call cg_error_exit_f

     ii = 0
     ! Write out the zones:
     do i=1,size(zones)

        sizes(:) = 0
        sizes(1) = zones(i)%il
        sizes(2) = zones(i)%jl
        sizes(3) = zones(i)%kl
        sizes(4) = zones(i)%il-1
        sizes(5) = zones(i)%jl-1
        sizes(6) = zones(i)%kl-1
        curSize = sizes(1)*sizes(2)*sizes(3)
        allocate(coorX(curSize), coorY(curSize), coorZ(curSize))

        ! Extract the coordinates
        do jj=1,curSize
           coorX(jj) = xx(3*ii + 3*jj-2)
           coorY(jj) = xx(3*ii + 3*jj-1)
           coorZ(jj) = xx(3*ii + 3*jj  )
        end do
     
        ii = ii + curSize
        call cg_zone_write_f(cg, base, zones(i)%name, sizes, Structured, zoneID, ierr)
        if (ierr .eq. CG_ERROR) call cg_error_exit_f

        ! Write the grid coordinates
        call cg_coord_write_f(cg, base, i, realDouble, 'CoordinateX', coorX, coordID, ierr)
        if (ierr .eq. CG_ERROR) call cg_error_exit_f

        call cg_coord_write_f(cg, base, i, realDouble, 'CoordinateY', coorY, coordID, ierr)
        if (ierr .eq. CG_ERROR) call cg_error_exit_f

        call cg_coord_write_f(cg, base, i, realDouble, 'CoordinateZ', coorZ, coordID, ierr)
        if (ierr .eq. CG_ERROR) call cg_error_exit_f
        
        deallocate(coorX, coorY, coorZ)
        
        ! Write any boundary conditions
        do boco=1,size(zones(i)%bocos)
           call cg_boco_write_f(cg, base, i, trim(zones(i)%bocos(boco)%name), &
                zones(i)%bocos(boco)%type, PointRange, &
                2, zones(i)%bocos(boco)%ptRange, BCout, ierr)
           if (ierr .eq. CG_ERROR) call cg_error_exit_f

           call cg_goto_f(cg, base, ierr, 'Zone_t', i, "ZoneBC_t", 1, &
                "BC_t", BCOut, "end")
           if (ierr .eq. CG_ERROR) call cg_error_exit_f

           if (zones(i)%bocos(boco)%family .ne. "") then 
              call cg_famname_write_f(zones(i)%bocos(boco)%family, ierr)
              if (ierr .eq. CG_ERROR) call cg_error_exit_f
           end if
        end do

        ! Write any b2b conditions
        do iB2B=1,size(zones(i)%b2bs)
           call cg_1to1_write_f(cg, base, i, &
                zones(i)%b2bs(iB2b)%name, &
                zones(i)%b2bs(iB2b)%donorName, &
                zones(i)%b2bs(iB2b)%ptRange, &
                zones(i)%b2bs(iB2b)%donorRange, &
                zones(i)%b2bs(iB2b)%transform, nCon, ierr)
           if (ierr .eq. CG_ERROR) call cg_error_exit_f
        end do
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
end subroutine writeStructuredCGNS

