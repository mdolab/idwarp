! ====================================================================
! File: allocateGridDoms.f90
! Author: C.A.(Sandy) Mader
! Date Started: July 29, 2014
! Date Modified:

subroutine allocateGridDoms(nZonesIn)
  
  use precision
  use gridData
  implicit none

  ! Subroutine Variables
  integer(kind=intType):: nZonesIn
  
  ! local Variables
  integer(kind=intType)::ierr

  ! begin execution
  
  if(.not. allocated(gridDoms))then
     allocate(gridDoms(nZonesIn),STAT=ierr)
  end if
  nZones = nZonesIn
end subroutine allocateGridDoms
