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
  !print *,'computing Nodal Properties',initialPoint,hasSymmetry,symDir

  ! First loop over the Farfield boundaries
  do pt = 1, nUniqueBoundaryPoints
     ! for each node loop over the connected elements and sum the area
     ! and area weighted normal
     zone = uniqueBoundaryNodes(pt)%zone
     sumArea = 0
     sumNormal = 0
     !surfSecCounter = uniqueBoundaryNodes(pt)%surfSection
     do elem = 1,30 ! This is hardcoded for the moment to reduce allocate calls
        elemIdx = uniqueBoundaryNodes(pt)%connectedElements(elem,1)
        surfSecCounter = uniqueBoundaryNodes(pt)%connectedElements(elem,2)
        if (elemIdx .ge. 0) then
           !print *,'elemIdx',elem,elemIdx,surfSecCounter,shape( gridDoms(zone)%surfaceSections(surfSecCounter)%elemAreaMag),shape(gridDoms(zone)%surfaceSections(surfSecCounter)%elemNormal)
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
     if(hasSymmetry)then
        uniqueBoundaryNodes(pt)%symAi = sumArea

        if(initialPoint)then
           uniqueBoundaryNodes(pt)%symNormal0 = sumNormal!*(-1*symDir)
           uniqueBoundaryNodes(pt)%symNormal0(iSymm) = uniqueBoundaryNodes(pt)%symNormal0(iSymm)*(-1)
           uniqueBoundaryNodes(pt)%symloc0 = uniqueBoundaryNodes(pt)%loc0!*(-1*symDir)
           uniqueBoundaryNodes(pt)%symloc0(iSymm) = uniqueBoundaryNodes(pt)%symloc0(iSymm)*-1
        else
           uniqueBoundaryNodes(pt)%symNormal = sumNormal!*(-1*symDir)
           uniqueBoundaryNodes(pt)%symNormal(iSymm) = uniqueBoundaryNodes(pt)%symNormal(iSymm)*(-1) 
           uniqueBoundaryNodes(pt)%symloc = uniqueBoundaryNodes(pt)%loc!*(-1*symDir)
           uniqueBoundaryNodes(pt)%symloc(iSymm) =  uniqueBoundaryNodes(pt)%symloc(iSymm)*(-1)
           v1 = uniqueBoundaryNodes(pt)%symNormal0
           v2 = uniqueBoundaryNodes(pt)%symNormal
           call getRotationMatrix3d(v1,v2,Mi)
           uniqueBoundaryNodes(pt)%symMi = Mi
           bi = (uniqueBoundaryNodes(pt)%symloc-uniqueBoundaryNodes(pt)%symloc0)*(-1*symDir)
           uniqueBoundaryNodes(pt)%symbi = bi
        end if
     end if
  end do


  ! Then loop over the surface boundaries
  do pt = 1, nUniqueSurfPoints
     !print *,'in pt',pt
     ! for each node loop over the connected elements and sum the area
     ! and area weighted normal
     zone = uniqueSurfaceNodes(pt)%zone
     sumArea = 0
     sumNormal = 0
     !surfSecCounter = uniqueSurfaceNodes(pt)%surfSection
     do elem = 1,30 ! This is hardcoded for the moment to reduce allocate calls
        elemIdx = uniqueSurfaceNodes(pt)%connectedElements(elem,1)
        surfSecCounter = uniqueSurfaceNodes(pt)%connectedElements(elem,2)
        if (elemIdx .ge. 0) then
           ! this is a connected element
           area = gridDoms(zone)%surfaceSections(surfSecCounter)%elemAreaMag(elemIdx)
           normal = gridDoms(zone)%surfaceSections(surfSecCounter)%elemNormal(elemIdx,:)
           !print *,'sumArea',elem,sumArea,area
           sumArea = sumArea + area
           !print *,'sumNormal',elem,sumNormal,normal
           sumNormal = sumNormal + area*normal
        end if
     end do
     sumNormal = sumNormal/sumArea
     uniqueSurfaceNodes(pt)%Ai = sumArea
     if(initialPoint)then
        uniqueSurfaceNodes(pt)%normal0 = sumNormal
     else
        !print *,'pt',pt
        uniqueSurfaceNodes(pt)%normal = sumNormal
        v1 = uniqueSurfaceNodes(pt)%normal0
        v2 = uniqueSurfaceNodes(pt)%normal
        call getRotationMatrix3d(v1,v2,Mi)
        uniqueSurfaceNodes(pt)%Mi = Mi
        bi = uniqueSurfaceNodes(pt)%loc-uniqueSurfaceNodes(pt)%loc0
        uniqueSurfaceNodes(pt)%bi = bi

     end if
     if(hasSymmetry)then
        uniqueSurfaceNodes(pt)%symAi = sumArea

        if(initialPoint)then
           uniqueSurfaceNodes(pt)%symNormal0 = sumNormal!*(-1*symDir)
           uniqueSurfaceNodes(pt)%symNormal0(iSymm) =uniqueSurfaceNodes(pt)%symNormal0(iSymm)*(-1)
           uniqueSurfaceNodes(pt)%symLoc0 = uniqueSurfaceNodes(pt)%loc0!*(-1*symDir)
           uniqueSurfaceNodes(pt)%symLoc0(iSymm) = uniqueSurfaceNodes(pt)%symLoc0(iSymm)*(-1)
           ! print *,'symloc',uniqueSurfaceNodes(pt)%symLoc0,uniqueSurfaceNodes(pt)%loc0,(-1*symDir)
           ! stop
        else
           uniqueSurfaceNodes(pt)%symNormal = sumNormal!*(-1*symDir)
           uniqueSurfaceNodes(pt)%symNormal(iSymm) = uniqueSurfaceNodes(pt)%symNormal(iSymm)*(-1) 
           uniqueSurfaceNodes(pt)%symLoc = uniqueSurfaceNodes(pt)%loc!*(-1*symDir)
           uniqueSurfaceNodes(pt)%symLoc(iSymm) = uniqueSurfaceNodes(pt)%symLoc(iSymm)*(-1)
           v1 = uniqueSurfaceNodes(pt)%symNormal0
           v2 = uniqueSurfaceNodes(pt)%symNormal
           call getRotationMatrix3d(v1,v2,Mi)
           uniqueSurfaceNodes(pt)%symMi = Mi
           bi = (uniqueSurfaceNodes(pt)%symloc-uniqueSurfaceNodes(pt)%symloc0)!*(-1*symDir)
           uniqueSurfaceNodes(pt)%symbi = bi
        end if
     end if

  end do

  
end subroutine computeNodalProperties
