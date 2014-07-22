! ====================================================================
! File: computeWi.f90
! Author: C.A.(Sandy) Mader
! Date Started: July 15, 2014
! Date Modified:

subroutine computeSi(surfZone,surfSec,surfNode,r,si)

  use precision
  use gridData
  implicit none

  ! Subroutine Variables
  integer(kind=intType)::surfZone,surfSec,surfNode
  real(kind=realType), dimension(physDim)::r,si

  ! local variables
  real(kind=realType), dimension(physDim)::bi
  real(kind=realType), dimension(physDim,physDim)::Mi

  Mi = gridDoms(surfZone)%surfaceSections(surfSec)%nodalRotation(surfNode)
  bi =  gridDoms(surfZone)%surfaceSections(surfSec)%nodalDisplacement(surfNode)

  Si = matMul(Mi,r)+bi-r

end subroutine computeSi
