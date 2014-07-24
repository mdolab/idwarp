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
  integer(kind=intType):: zone,sec,node,pt
  integer(kind=intType):: volSecCounter, surfSecCounter,globalIndex
  real(kind=realType), dimension(physDim):: numerator
  real(kind=realType)::  denomenator
  real(kind=realType), dimension(physDim)::r,si
  real(kind=realType):: wi
  real(kind=realType):: timeA,timeB

  ! Begin Execution
  call cpu_time(timeA)
  aexp = 3.0
  bexp = 5.0
  Ldef = 15.0
  alpha = 0.1250

  !loop over the volume nodes
  do zone = 1,nZones
     volSecCounter = 1
     do sec = 1,gridDoms(zone)%nSections 
        if(gridDoms(zone)%isVolumeSection(sec))then
           do node = 1,gridDoms(zone)%nVertices
              !print *,'warping node',node,gridDoms(zone)%nVertices,sec,zone
              if(.not. gridDoms(zone)%isSurfaceNode(node))then
                 numerator = 0
                 denomenator = 0
                 r = gridDoms(zone)%points0(node,:)
                 ! boundary not needed because Si is always 0
                 ! ! First loop over the Farfield boundaries
                 ! do pt = 1, nUniqueBoundaryPoints
                 !    call computeWi(pt,.false.,r,wi)
                 !    call computeSi(pt,.false.,r,Si)
                 !    numerator = numerator + Wi*Si
                 !    denomenator = denomenator + Wi
                 ! end do
                 ! now loop over surfaces
                 do pt = 1, nUniqueSurfPoints
                    call computeWi(pt,.true.,r,wi)
                    call computeSi(pt,.true.,r,Si)
                    numerator = numerator + Wi*Si
                    denomenator = denomenator + Wi
                 end do
              end if
              gridDoms(zone)%dx(node,:) = numerator/denomenator
              gridDoms(zone)%points(node,:)=gridDoms(zone)%points0(node,:)+gridDoms(zone)%dx(node,:)
           end do
        end if
     end do
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
     gridDoms(zone)%points(globalIndex,:) = uniqueSurfaceNodes(pt)%loc  
  end do
end subroutine updateVolumeCoordinates
