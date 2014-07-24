! ====================================================================
! File: getSurfaceElementCenters.f90
! Author: C.A.(Sandy) Mader
! Date Started: July 14, 2014
! Date Modified:

subroutine getSurfaceElementCenterAndArea()

  use precision
  use gridData
  implicit none

  ! Local Variable
  integer(kind=intType)::zone,sec,ierr
  real(kind=realType),dimension(:,:), allocatable :: points
  real(kind=realType),dimension(:), allocatable:: center,area,normal
  real(kind=realType):: areaMag
  
  integer(kind=intType):: i, elem
  integer(kind=intType):: surfSecCounter,npts,nDim,pointIdx

  ! Begin Execution

  print *,' in get surface element centers',nZones
  ! loop over the zones and sections.
  do zone = 1,nZones
     surfSecCounter = 1
     print *,'nSections',gridDoms(zone)%nSections
     do sec = 1,gridDoms(zone)%nSections
        print *,'isSurface',sec,(.not. gridDoms(zone)%isVolumeSection(sec))
        if(.not. gridDoms(zone)%isVolumeSection(sec))then
           ! For each surface section, compute the element center and store
           npts = gridDoms(zone)%surfaceSections(surfSecCounter)%nConn
           nDim = physDim
           allocate(points(npts,nDim),center(nDim),area(nDim),STAT=ierr)
           allocate(normal(nDim),STAT=ierr)
           do elem =1,gridDoms(zone)%surfaceSections(surfSecCounter)%nElem
              do i = 1,nPts
                 pointIdx = gridDoms(zone)%surfaceSections(surfSecCounter)%elements(i,elem)
                 points(i,:) = gridDoms(zone)%points(pointIdx,:)
              end do
              call getElementCenter(nPts,nDim,points,center)
              call getElementArea(nPts,nDim,points,center,area,areaMag,normal)

              gridDoms(zone)%surfaceSections(surfSecCounter)%elemCenter(elem,:) = center
              gridDoms(zone)%surfaceSections(surfSecCounter)%elemArea(elem,:) = area
              gridDoms(zone)%surfaceSections(surfSecCounter)%elemAreaMag(elem) = areaMag
              gridDoms(zone)%surfaceSections(surfSecCounter)%elemNormal(elem,:) = normal

           end do
           deallocate(points,center,area,normal,STAT=ierr)
           surfSecCounter = surfSecCounter + 1
        end if
     end do
  end do

end subroutine getSurfaceElementCenterAndArea
