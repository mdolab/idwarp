subroutine readCGNS(cgns_file)
  use gridData
  use gridInput
  use communication
  use CGNSGrid

  implicit none

  ! Input Arguments
  character*(*),intent(in) :: cgns_file

  ! Working
  integer(kind=intType) :: cg, ierr, i
  integer(kind=intType):: CellDim, PhysDim, nZones, base, nbases
  integer(kind=intType) :: nstructured, nunstructured, zoneType
  character(len=32) :: baseName

  ! Set the default family names
  defaultFamName(BCAxisymmetricWedge) = 'axi'
  defaultFamName(BCDegenerateLine) = 'degenerate'
  defaultFamName(BCDegeneratePoint) ='degenerate'
  defaultFamName(BCDirichlet) = 'dirichlet'
  defaultFamName(BCExtrapolate) = 'extrap'
  defaultFamName(BCFarfield) = 'far'
  defaultFamName(BCGeneral) = 'general'
  defaultFamName(BCInflow) = 'inflow'
  defaultFamName(BCInflowSubsonic) = 'inflow'
  defaultFamName(BCInflowSupersonic) = 'inflow'
  defaultFamName(BCNeumann) = 'neumann'
  defaultFamName(BCOutflow) = 'outflow'
  defaultFamName(BCOutflowSubsonic) = 'outflow'
  defaultFamName(BCOutflowSupersonic)  ='outflow'
  defaultFamName(BCSymmetryPlane) = 'sym'
  defaultFamName(BCSymmetryPolar) = 'sympolar'
  defaultFamName(BCTunnelInflow) = 'inflow'
  defaultFamName(BCTunnelOutflow) = 'outflow'
  defaultFamName(BCWall) = 'wall'
  defaultFamName(BCWallInviscid) = 'wall'
  defaultFamName(BCWallViscous) = 'wall'
  defaultFamName(BCWallViscousHeatFlux) = 'wall'
  defaultFamName(BCWallViscousIsothermal) = 'wall'
  defaultFamName(FamilySpecified) = 'wall'

  ! Do the I/O that is common to both types of grids

  if (myid == 0) then
     print *, ' -> Reading CGNS File: ', cgns_file

     ! Open and get the number of zones:
     call cg_open_f(trim(cgns_file), CG_MODE_READ, cg, ierr)
     if (ierr .eq. CG_ERROR) call cg_error_exit_f

     call cg_nbases_f(cg, nbases, ierr)
     if (ierr .eq. CG_ERROR) call cg_error_exit_f

     if (nbases .gt. 1) then
        print *, ' ** Warning: IDWarp only reads the first base in a cgns file'
     end if

     base = 1_intType

     call cg_base_read_f(cg, base, basename, CellDim, PhysDim, ierr)
     if (ierr .eq. CG_ERROR) call cg_error_exit_f
     if (cellDim .ne. 3 .or. PhysDim .ne. 3) then
        print *, 'The Cells must 3 dimensional'
        stop
     end if

     call cg_nzones_f(cg, base, nZones, ierr);
     if (ierr .eq. CG_ERROR) call cg_error_exit_f
     print *, '   -> Number of Zones:', nzones

     ! Determine if we have structured or unstructured zones. We can
     ! only deal with one or the other.
     nStructured = 0
     nUnstructured = 0
     do i=1, nZones
        call cg_zone_type_f(cg, base, i, zoneType, ierr)
        if (zoneType == Structured) then
           nStructured = nStructured + 1
        else if (zoneType == Unstructured) then
           nUnstructured = nUnstructured + 1
        end if
     end do
  end if

  ! Broadcast the results to everyone
  call MPI_bcast(nStructured, 1, MPI_INTEGER, 0, warp_comm_world, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call MPI_bcast(nUnstructured, 1, MPI_INTEGER, 0, warp_comm_world, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! There are 4 possibilities:
  if (nStructured /= 0 .and. nUnstructured == 0) then
     cgnsStructured = .True.
     call readStructuredCGNS(cg)
  else if (nStructured == 0 .and. nUnstructured /= 0) then
     cgnsStructured = .False.
     call readUnstructuredCGNS(cg)
  else
     print *, "Error reading CGNS file. There are either no zones *OR* &
          &there are mixed structured and unstructured zones that &
          &cannot currently be handled."
     stop
  end if
end subroutine readCGNS
