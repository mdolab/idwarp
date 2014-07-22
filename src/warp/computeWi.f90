! ====================================================================
! File: computeWi.f90
! Author: C.A.(Sandy) Mader
! Date Started: July 15, 2014
! Date Modified:

subroutine computeWi(surfZone,surfSec,surfNode,r,wi)

  use precision
  use gridData
  implicit none

  ! Subroutine Variables
  integer(kind=intType)::surfZone,surfSec,surfNode
  real(kind=realType), dimension(physDim)::r
  real(kind=realType)::wi

  ! local Variables
  real(kind=realType), dimension(physDim)::ri
  real(kind=realType):: Ai

  
  ! Begin Execution

  Ai = gridDoms(surfZone)%surfaceSections(surfSec)%nodalArea(surfNode)
  ri =  gridDoms(surfZone)%points(surfNode,:)
  call getDistance(r,ri,dist)
  wi = Ai *((Ldef/dist)**aExp+(alpha*Ldef/dist)**bExp)

end subroutine computeWi
