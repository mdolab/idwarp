! ====================================================================
! File: computeWi.f90
! Author: C.A.(Sandy) Mader
! Date Started: July 15, 2014
! Date Modified:

subroutine computeSi(surfIdx,isSurf,r,si)

  use precision
  use gridData
  implicit none

  ! Subroutine Variables
  integer(kind=intType)::surfIdx
  real(kind=realType), dimension(3)::r,si
  logical :: isSurf

  ! local variables
  !real(kind=realType), dimension(3)::v1,v2
  real(kind=realType), dimension(3)::bi
  real(kind=realType), dimension(3,3)::Mi
  
  ! if(isSurf)then
  !    v1 = uniqueSurfaceNodes(surfIdx)%normal0
  !    v2 = uniqueSurfaceNodes(surfIdx)%normal
  !    bi =  uniqueSurfaceNodes(surfIdx)%loc-uniqueSurfaceNodes(surfIdx)%loc0
  ! else
  !    v1 = uniqueBoundaryNodes(surfIdx)%normal0
  !    v2 = uniqueBoundaryNodes(surfIdx)%normal
  !    bi = uniqueBoundaryNodes(surfIdx)%loc-uniqueBoundaryNodes(surfIdx)%loc0
  ! end if

  ! call getRotationMatrix3d(v1,v2,Mi)
! if(isSurf)then
     bi =  uniqueSurfaceNodes(surfIdx)%bi
     Mi =  uniqueSurfaceNodes(surfIdx)%Mi
!  print *,'bi,mi',bi,mi
  ! else
  !    bi = uniqueBoundaryNodes(surfIdx)%bi
  !    Mi = uniqueBoundaryNodes(surfIdx)%Mi
  ! end if
  
  Si = matMul(Mi,r)+bi-r
  !print *,'Si',Si

end subroutine computeSi

subroutine computeSiSymm(surfIdx,isSurf,r,si)

  use precision
  use gridData
  implicit none

  ! Subroutine Variables
  integer(kind=intType)::surfIdx
  real(kind=realType), dimension(3)::r,si
  logical :: isSurf

  ! local variables
  !real(kind=realType), dimension(3)::v1,v2
  real(kind=realType), dimension(3)::bi
  real(kind=realType), dimension(3,3)::Mi
  
  ! call getRotationMatrix3d(v1,v2,Mi)
!  if(isSurf)then
     bi =  uniqueSurfaceNodes(surfIdx)%symbi
     Mi =  uniqueSurfaceNodes(surfIdx)%symMi
  ! else
  !    bi = uniqueBoundaryNodes(surfIdx)%symbi
  !    Mi = uniqueBoundaryNodes(surfIdx)%symMi
  ! end if
  
  Si = matMul(Mi,r)+bi-r

end subroutine computeSiSymm
