! ====================================================================
! File: setGetVolumeCoordinates.f90
! Author: C.A.(Sandy) Mader
! Date Started: July 29, 2014
! Date Modified:

subroutine setVolumeCoordinates(zone,nZonesIn,nVertices,vertices,initCall)

  use precision
  use gridData
  implicit none

  ! Subroutine Variables
  integer(kind=intType):: zone,nZonesIn,nVertices
  real(kind=realType),dimension(nVertices,3)::vertices
  logical:: initCall

  ! Local Variables
  integer(kind=intType):: ierr, i
  physDim = 3
  ! Begin Execution
  if(initCall)then

     if(.not. allocated(gridDoms)) then 
        call allocateGridDoms(nZonesIn)
     end if

     if(.not. allocated(gridDoms(zone)%points))then

        allocate(gridDoms(zone)%points(nVertices,3),STAT=ierr)

        allocate(gridDoms(zone)%points0(nVertices,3),STAT=ierr)

        allocate(gridDoms(zone)%dx(nVertices,3),STAT=ierr)

        allocate(gridDoms(zone)%isSurfaceNode(nVertices),STAT=ierr)

     end if

     gridDoms(zone)%isSurfaceNode=.False.
     gridDoms(zone)%points0 = vertices

  end if
  ! do i = 1,nVertices
  !    print *,'x',i,vertices(i,:)
  ! end do
  gridDoms(zone)%nVertices = nVertices
  gridDoms(zone)%points = vertices
  gridDoms(zone)%dx = 0

end subroutine setVolumeCoordinates

subroutine getVolumeCoordinates(zone,nVertices,vertices)

  use precision
  use gridData
  implicit none

  ! Subroutine Variables
  integer(kind=intType)::zone,nVertices
  real(kind=realType),dimension(nVertices,3),intent(out)::vertices
  
  ! begin Execution
  vertices = gridDoms(zone)%points

end subroutine getVolumeCoordinates

subroutine getNVertices(zone,nVertices)
  use precision
  use gridData
  implicit none

  ! Subroutine Variables
  integer(kind=intType),intent(in)::zone
  integer(kind=intType),intent(out)::nVertices

  ! begin execution

  nVertices = gridDoms(zone)%nVertices

end subroutine getNVertices
  
