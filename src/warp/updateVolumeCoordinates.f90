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



  ! Begin Execution

  aexp = 3.0
  bexp = 5.0
  L = 15.0
  alpha = 0.1250

  !loop over the volume nodes
  do zone = 1,nZones
     volSecCounter = 1
     do sec = 1,gridDoms(zone)%nSections 
        if(gridDoms(zone)%isVolumeSection(sec))then
           do node = 1,gridDoms(zone)%nVertices
              if(.not. gridDoms(zone)%isSurfaceNode(node))then
                 numerator = 0
                 denomenator = 0
                 r = gridDoms(zone)%points(node,:)
                 do surfZone = 1,nZones
                    surfSecCounter = 1
                    do surfSec = 1,gridDoms(surfZone)%nSections 
                       if(.not. gridDoms(surfZone)%isVolumeSection(surfSec))then
                          do surfNode = 1,gridDoms(surfZone)%nSurfNodes
                             call computeWi(surfZone,surfSec,surfNode,r,wi)
                             call computeSi(surfZone,surfSec,surfNode,r,Si)
                             numerator = numerator + Wi*Si
                             denomenator = denomenator + Wi
                          end do
                       end if
                    end do
                 end do
                 gridDoms(zone)%dx(node,:) = numerator/denomenator
              end if
           end do
        end if
     end do
  end do

end subroutine updateVolumeCoordinates
