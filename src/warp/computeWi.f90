! ====================================================================
! File: computeWi.f90
! Author: C.A.(Sandy) Mader
! Date Started: July 15, 2014
! Date Modified:

subroutine computeWi(surfIdx,isSurf,r,wi)

  use precision
  use gridData
  implicit none

  ! Subroutine Variables
  integer(kind=intType)::surfIdx
  real(kind=realType), dimension(3)::r
  real(kind=realType)::wi
  logical :: isSurf

  ! local Variables
  real(kind=realType), dimension(3)::ri
  real(kind=realType):: Ai,dist

  
  ! Begin Execution
  !print *,'in wi',surfIdx,isSurf,r

!  if(isSurf)then

     Ai = uniqueSurfaceNodes(surfIdx)%Ai
     ri = uniqueSurfaceNodes(surfIdx)%loc
     !print *,'Ai', Ai,uniqueSurfaceNodes(surfIdx)%Ai
  ! else
  !    Ai = uniqueBoundaryNodes(surfIdx)%Ai
  !    ri = uniqueBoundaryNodes(surfIdx)%loc
  !    !print *,'ai b',Ai,uniqueBoundaryNodes(surfIdx)%Ai
  ! end if
!  print *,'ri',ri
  call getDistance(r,ri,dist)
  ! print *,'dist',dist
  ! print *,'Ai',Ai
  ! print *,'ldef',ldef
  ! print *,'aexp',aexp
  ! print *,'bexp',bexp
  ! print *,'alpha',alpha
  ! stop
!  print *,'wi',Ai,lDef,dist,aexp,'t2',alpha,bexp
  wi = Ai *((Ldef/dist)**aExp+(alpha*Ldef/dist)**bExp)

end subroutine computeWi

subroutine computeWiSymm(surfIdx,isSurf,r,wi)

  use precision
  use gridData
  implicit none

  ! Subroutine Variables
  integer(kind=intType)::surfIdx
  real(kind=realType), dimension(3)::r
  real(kind=realType)::wi
  logical :: isSurf

  ! local Variables
  real(kind=realType), dimension(3)::ri
  real(kind=realType):: Ai,dist

  
  ! Begin Execution
  !print *,'in wi',surfIdx,isSurf,r

  ! if(isSurf)then

  Ai = uniqueSurfaceNodes(surfIdx)%Ai
  ri = uniqueSurfaceNodes(surfIdx)%symloc
  !    !print *,'Ai', Ai,uniqueSurfaceNodes(surfIdx)%Ai
  ! else
  !Ai = uniqueBoundaryNodes(surfIdx)%Ai
  !ri = uniqueBoundaryNodes(surfIdx)%symloc
  !    !print *,'ai b',Ai,uniqueBoundaryNodes(surfIdx)%Ai
  ! end if
!  print *,'ri',ri
  call getDistance(r,ri,dist)
  ! print *,'dist',dist
  ! print *,'Ai',Ai
  ! print *,'ldef',ldef
  ! print *,'aexp',aexp
  ! print *,'bexp',bexp
  ! print *,'alpha',alpha
  ! stop
!  print *,'wi',Ai,lDef,dist,aexp,'t2',alpha,bexp
  wi = Ai *((Ldef/dist)**aExp+(alpha*Ldef/dist)**bExp)

end subroutine computeWiSymm
