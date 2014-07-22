! ====================================================================
! File: computeNodalProperties.f90
! Author: C.A.(Sandy) Mader
! Date Started: July 14, 2014
! Date Modified:

subroutine computeNodalProperties()

  use precision
  use gridData
  implicit none

  ! Local Variable
  integer(kind=intType)::zone,sec,ierr





  ! begin Execution

  !Loop over the surface nodes and compute their areas and normals
  do zone = 1,nZones
     surfSecCounter = 1 
     do sec = 1,gridDoms(zone)%nSections
        print *,'nodal properties: section',sec
        if(.not. gridDoms(zone)%isVolumeSection(sec))then
           nPts = gridDoms(zone)%surfaceSections(surfSecCounter)%nSurf

           do node = 1,nPts
              ! for each node loop over the connected elements and sum the area
              ! and area weighted normal
              sumArea = 0.
              sumNormal = 0.
              do elem = 1,10 ! This is hardcoded for the moment to reduce allocate calls
                 elemIdx = gridDoms(zone)%surfaceSections(surfSecCounter)%uniqueSurfaceNodes(node)%connectedElements(elem)
                 if (elemIdx .ge. 0) then
                    ! this is a connected element
                    area = gridDoms(zone)%surfaceSections(surfSecCounter)%elemAreaMag(elemIdx)
                    normal = gridDoms(zone)%surfaceSections(surfSecCounter)%normal(elemIdx,:)
                    sumArea = sumArea + area
                    sumNormal = sumNormal + area*normal
                 end if
              end do
              sumNormal = sumNormal/sumArea
              gridDoms(zone)%surfaceSections(surfSecCounter)%uniqueSurfaceNodes(node)%Ai = sumArea
              gridDoms(zone)%surfaceSections(surfSecCounter)%uniqueSurfaceNodes(node)%Normal = sumNormal
           end do
        end if
     end do
  end do

end subroutine computeNodalProperties
