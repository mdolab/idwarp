! ====================================================================
! File: computeNodalProperties.f90
! Author: C.A.(Sandy) Mader
! Date Started: July 14, 2014
! Date Modified:

subroutine computeNodalProperties(initialPoint)

  use precision
  use gridData
  implicit none

  ! Subroutine Variable
  logical :: initialPoint

  ! Local Variable
  integer(kind=intType)::zone,sec,ierr,elem,elemIdx,pt,surfSecCounter
  real(kind=realType):: area,sumArea
  real(kind=realType), dimension(3)::normal,sumNormal,v1,v2, bi
  real(kind=realType), dimension(3,3)::Mi
  ! begin Execution
  print *,'computing Nodal Properties',initialPoint
  ! First loop over the Farfield boundaries
  do pt = 1, nUniqueBoundaryPoints
     ! for each node loop over the connected elements and sum the area
     ! and area weighted normal
     zone = uniqueBoundaryNodes(pt)%zone
     sumArea = 0
     sumNormal = 0
     !surfSecCounter = uniqueBoundaryNodes(pt)%surfSection
     do elem = 1,10 ! This is hardcoded for the moment to reduce allocate calls
        elemIdx = uniqueBoundaryNodes(pt)%connectedElements(elem,1)
        surfSecCounter = uniqueBoundaryNodes(pt)%connectedElements(elem,2)
        if (elemIdx .ge. 0) then
           ! this is a connected element
           area = gridDoms(zone)%surfaceSections(surfSecCounter)%elemAreaMag(elemIdx)
           normal = gridDoms(zone)%surfaceSections(surfSecCounter)%elemNormal(elemIdx,:)
           sumArea = sumArea + area
           sumNormal = sumNormal + area*normal
        end if
     end do
     sumNormal = sumNormal/sumArea
     uniqueBoundaryNodes(pt)%Ai = sumArea

     if(initialPoint)then
        uniqueBoundaryNodes(pt)%normal0 = sumNormal
     else
        uniqueBoundaryNodes(pt)%normal = sumNormal
        v1 = uniqueBoundaryNodes(pt)%normal0
        v2 = uniqueBoundaryNodes(pt)%normal
        call getRotationMatrix3d(v1,v2,Mi)
        uniqueBoundaryNodes(pt)%Mi = Mi
        bi = uniqueBoundaryNodes(pt)%loc-uniqueBoundaryNodes(pt)%loc0
        uniqueBoundaryNodes(pt)%bi = bi
     end if

  end do


  ! Then loop over the surface boundaries
  do pt = 1, nUniqueSurfPoints
     ! for each node loop over the connected elements and sum the area
     ! and area weighted normal
     zone = uniqueSurfaceNodes(pt)%zone
     sumArea = 0
     sumNormal = 0
     !surfSecCounter = uniqueSurfaceNodes(pt)%surfSection
     do elem = 1,10 ! This is hardcoded for the moment to reduce allocate calls
        elemIdx = uniqueSurfaceNodes(pt)%connectedElements(elem,1)
        surfSecCounter = uniqueSurfaceNodes(pt)%connectedElements(elem,2)
        if (elemIdx .ge. 0) then
           ! this is a connected element
           !print *,'areaMag',zone,surfSecCounter,elemIdx,shape(gridDoms(zone)%surfaceSections(surfSecCounter)%elemAreaMag)
           area = gridDoms(zone)%surfaceSections(surfSecCounter)%elemAreaMag(elemIdx)
           normal = gridDoms(zone)%surfaceSections(surfSecCounter)%elemNormal(elemIdx,:)
           sumArea = sumArea + area
           sumNormal = sumNormal + area*normal
        end if
     end do
     sumNormal = sumNormal/sumArea
     uniqueSurfaceNodes(pt)%Ai = sumArea
     if(initialPoint)then
        uniqueSurfaceNodes(pt)%normal0 = sumNormal
     else
        uniqueSurfaceNodes(pt)%normal = sumNormal
        v1 = uniqueSurfaceNodes(pt)%normal0
        v2 = uniqueSurfaceNodes(pt)%normal
        call getRotationMatrix3d(v1,v2,Mi)
        uniqueSurfaceNodes(pt)%Mi = Mi
        bi = uniqueSurfaceNodes(pt)%loc-uniqueSurfaceNodes(pt)%loc0
        uniqueSurfaceNodes(pt)%bi = bi
     end if

  end do

  ! !Loop over the surface nodes and compute their areas and normals
  ! do zone = 1,nZones
  !    surfSecCounter = 1 
  !    do sec = 1,gridDoms(zone)%nSections
  !       print *,'nodal properties: section',sec
  !       if(.not. gridDoms(zone)%isVolumeSection(sec))then
  !          nPts = gridDoms(zone)%surfaceSections(surfSecCounter)%nSurf

  !          do node = 1,nPts
  !             ! for each node loop over the connected elements and sum the area
  !             ! and area weighted normal
  !             sumArea = 0.
  !             sumNormal = 0.
  !             do elem = 1,10 ! This is hardcoded for the moment to reduce allocate calls
  !                elemIdx = gridDoms(zone)%surfaceSections(surfSecCounter)%uniqueSurfaceNodes(node)%connectedElements(elem)
  !                if (elemIdx .ge. 0) then
  !                   ! this is a connected element
  !                   area = gridDoms(zone)%surfaceSections(surfSecCounter)%elemAreaMag(elemIdx)
  !                   normal = gridDoms(zone)%surfaceSections(surfSecCounter)%normal(elemIdx,:)
  !                   sumArea = sumArea + area
  !                   sumNormal = sumNormal + area*normal
  !                end if
  !             end do
  !             sumNormal = sumNormal/sumArea
  !             gridDoms(zone)%surfaceSections(surfSecCounter)%uniqueSurfaceNodes(node)%Ai = sumArea
  !             gridDoms(zone)%surfaceSections(surfSecCounter)%uniqueSurfaceNodes(node)%Normal = sumNormal
  !          end do
  !       end if
  !    end do
  ! end do

end subroutine computeNodalProperties
