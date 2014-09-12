! ====================================================================
! File: updateVolumeCoordinates.f90
! Author: C.A.(Sandy) Mader
! Date Started: July 15, 2014
! Date Modified:

subroutine updateVolumeCoordinates()

  use precision
  use gridData
  implicit none


  ! Local Variables
  integer(kind=intType):: zone,sec,node,pt,nodeNumber,i
  integer(kind=intType):: volSecCounter, surfSecCounter,globalIndex
  real(kind=realType), dimension(3):: numerator
  real(kind=realType)::  denomenator
  real(kind=realType), dimension(3)::r,si
  real(kind=realType):: wi
  real(kind=realType):: timeA,timeB
  real(kind=realType), dimension(3)::bi
  real(kind=realType), dimension(3,3)::Mi
  real(kind=realType), dimension(3)::ri,dx
  real(kind=realType):: Ai,dist

  ! Begin Execution
  call cpu_time(timeA)
  aexp = 3.0
  bexp = 5.0
  Ldef = 3.5!15.0
  alpha = 1.0!0.1250

  nodeNumber = 10048!12944
  !loop over the volume nodes
  do zone = 1,nZones
     volSecCounter = 1
     ! do sec = 1,gridDoms(zone)%nSections 
     !    if(gridDoms(zone)%isVolumeSection(sec))then
     !print *,'in Update Volume',gridDoms(zone)%nVertices
           do node = 1,gridDoms(zone)%nVertices
              !print *,'warping node',node,gridDoms(zone)%nVertices,sec,zone
              !if(node .eq. nodeNumber)then
              if(.not. gridDoms(zone)%isSurfaceNode(node))then

                 numerator = 0
                 denomenator = 0
                 r = gridDoms(zone)%points0(node,:)
                 ! boundary Si not needed because it is always 0
                 ! just compute Wi
                 ! First loop over the Farfield boundaries
                 do pt = 1, nUniqueBoundaryPoints
                    !print *,'boundarypoint',pt
                    !call computeWi(pt,.false.,r,wi)
                    Ai = uniqueBoundaryNodes(pt)%Ai
                    ri = uniqueBoundaryNodes(pt)%loc0
                    dx = r-ri
                    dist = sqrt(dx(1)*dx(1)+dx(2)*dx(2)+dx(3)*dx(3))

                    wi = Ai *((Ldef/dist)**aExp+(alpha*Ldef/dist)**bExp)

                    !call computeSi(pt,.false.,r,Si)
                    !numerator = numerator + Wi*Si
                    denomenator = denomenator + Wi
                 end do
                 ! now loop over surfaces
                 !print *,'r',r
                 do pt = 1, nUniqueSurfPoints
                    !print *,'SurfacePoint',pt
                    !call computeWi(pt,.true.,r,wi)
                    Ai = uniqueSurfaceNodes(pt)%Ai
                    ri = uniqueSurfaceNodes(pt)%loc0
                    dx = r-ri
                    dist = sqrt(dx(1)*dx(1)+dx(2)*dx(2)+dx(3)*dx(3))
                    wi = Ai *((Ldef/dist)**aExp+(alpha*Ldef/dist)**bExp)

                    !call computeSi(pt,.true.,r,Si)
                    bi =  uniqueSurfaceNodes(pt)%bi
                    Mi =  uniqueSurfaceNodes(pt)%Mi
                    
                    do i=1,3
                       Si(i) = Mi(i,1)*r(1)+Mi(i,2)*r(2)+Mi(i,3)*r(3)+bi(i)-r(i)
                    end do
                          
                    !print *,pt,'Wi,Si',Wi,Si
                    numerator = numerator + Wi*Si
                    denomenator = denomenator + Wi
                    ! if (node .eq. 303) then
                    !    print *,'num',numerator,denomenator
                    !    print *,'weights',pt,wi,'s',si
                    ! end if
                 end do
                 if(hasSymmetry)then
                    do pt = 1, nUniqueBoundaryPoints
                       !print *,'symmBoundaryPoint',pt
                       !call computeWiSymm(pt,.false.,r,wi)
                       Ai = uniqueBoundaryNodes(pt)%Ai
                       ri = uniqueBoundaryNodes(pt)%symloc0
                       dx = r-ri
                       dist = sqrt(dx(1)*dx(1)+dx(2)*dx(2)+dx(3)*dx(3))
                       
                       wi = Ai *((Ldef/dist)**aExp+(alpha*Ldef/dist)**bExp)
                       !call computeSiSymm(pt,.false.,r,Si)
                       !numerator = numerator + Wi*Si
                       denomenator = denomenator + Wi
                    end do
                    do pt = 1, nUniqueSurfPoints
                       !print *,'symmSurfacePoint',pt
                       !call computeWiSymm(pt,.true.,r,wi)
                       Ai = uniqueSurfaceNodes(pt)%Ai
                       ri = uniqueSurfaceNodes(pt)%symloc0
                       dx = r-ri
                       dist = sqrt(dx(1)*dx(1)+dx(2)*dx(2)+dx(3)*dx(3))
                       wi = Ai *((Ldef/dist)**aExp+(alpha*Ldef/dist)**bExp)
                       
                       !call computeSiSymm(pt,.true.,r,Si)
                       bi =  uniqueSurfaceNodes(pt)%symbi
                       Mi =  uniqueSurfaceNodes(pt)%symMi
                    
                       do i=1,3
                          Si(i) = Mi(i,1)*r(1)+Mi(i,2)*r(2)+Mi(i,3)*r(3)+bi(i)-r(i)
                       end do
                       !print *,pt,'Wi,Si',Wi,Si
                       numerator = numerator + Wi*Si
                       denomenator = denomenator + Wi
                       ! if (node .eq. 303) then
                       !    print *,'numSym',numerator,denomenator
                       !    print *,'weightsSym',pt,wi,'s',si
                       ! end if
                    end do
                 end if
                 ! if (node .eq.303) then
                 !    print *,'num,dem,',node,numerator,denomenator,numerator/denomenator
                 !    stop
                 ! end if
                 ! if (numerator(2) .ne. 0) then
                 !    print *,'num,dem,',node,numerator,denomenator,numerator/denomenator
                 ! end if
                 ! stop
                 gridDoms(zone)%dx(node,:) = numerator/denomenator
              else
                 
                 gridDoms(zone)%dx(node,:) = 0.0
              end if
           !end if
              !print *,'dx', gridDoms(zone)%dx(node,:),'p0',gridDoms(zone)%points0(node,:),'p',gridDoms(zone)%points(node,:)
              ! if (gridDoms(zone)%dx(node,2) .ne. 0)then
              !    print *,'dx', gridDoms(zone)%dx(node,:),'p0',gridDoms(zone)%points0(node,:),'p',gridDoms(zone)%points(node,:)
              ! end if
              ! if (node .eq.nodeNumber) then
              !    print *,'node dx',node,gridDoms(zone)%dx(node,:)
              ! end if
              gridDoms(zone)%points(node,:)=gridDoms(zone)%points0(node,:)+gridDoms(zone)%dx(node,:)
           end do
     !    end if
     ! end do
  end do
  call cpu_time(timeB)
  print *,'MeshWarping Time',timeB-timeA
  ! now set the surface coordinates
  ! First loop over the Farfield boundaries
  do pt = 1, nUniqueBoundaryPoints
     zone = uniqueBoundaryNodes(pt)%zone
     globalIndex = uniqueBoundaryNodes(pt)%globalIndex
     gridDoms(zone)%points(globalIndex,:) = uniqueBoundaryNodes(pt)%loc
  end do
  ! Then loop over the surface boundaries
  do pt = 1, nUniqueSurfPoints
     zone = uniqueSurfaceNodes(pt)%zone
     globalIndex = uniqueSurfaceNodes(pt)%globalIndex
     !print *,'surf',pt,globalIndex,gridDoms(zone)%points(globalIndex,:),uniqueSurfaceNodes(pt)%loc
     gridDoms(zone)%points(globalIndex,:) = uniqueSurfaceNodes(pt)%loc  
  end do
end subroutine updateVolumeCoordinates
