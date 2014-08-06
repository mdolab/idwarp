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
  integer(kind=intType):: surfSecCounter,npts,nDim,pointIdx,elemIdx

  ! Begin Execution

  print *,' in get surface element centers',nZones
  ! loop over the zones and sections.
  do zone = 1,nZones
     surfSecCounter = 1
     print *,'nSections',gridDoms(zone)%nSections
     do sec = 1,gridDoms(zone)%nSections
        !print *,'isSurface',sec,(.not. gridDoms(zone)%isVolumeSection(sec))
        if(.not. gridDoms(zone)%isVolumeSection(sec))then
           ! For each surface section, compute the element center and store
           !npts = gridDoms(zone)%surfaceSections(surfSecCounter)%nConn
           nDim = physDim
           allocate(center(nDim),area(nDim),STAT=ierr)
           allocate(normal(nDim),STAT=ierr)
           do elem =1,gridDoms(zone)%surfaceSections(surfSecCounter)%nElem
              nPts = gridDoms(zone)%surfaceSections(surfSecCounter)%elemPtr(elem+1)-gridDoms(zone)%surfaceSections(surfSecCounter)%elemPtr(elem)
              allocate(points(npts,nDim),STAT=ierr)
              do i = 1,nPts
               
                 elemIdx = gridDoms(zone)%surfaceSections(surfSecCounter)%elemPtr(elem)
                 !print *,'elemidx',zone,elem,surfSecCounter,elemIdx+i,elemidx,i,shape( gridDoms(zone)%surfaceSections(surfSecCounter)%elemPtr),shape(gridDoms(zone)%surfaceSections(surfSecCounter)%elemConn)
                 pointIdx = gridDoms(zone)%surfaceSections(surfSecCounter)%elemConn(elemIdx+i)
                 points(i,:) = gridDoms(zone)%points(pointIdx,:)
              end do
              !print *,'elem',elem
              call getElementCenter(nPts,nDim,points,center)
              !print *,'points',points,'c',center
              call getElementArea(nPts,nDim,points,center,area,areaMag,normal)
              !print *,'area',points,'c',center,'a',area,'am',areaMag,'n',normal
              gridDoms(zone)%surfaceSections(surfSecCounter)%elemCenter(elem,:) = center
              gridDoms(zone)%surfaceSections(surfSecCounter)%elemArea(elem,:) = area
              gridDoms(zone)%surfaceSections(surfSecCounter)%elemAreaMag(elem) = areaMag
              gridDoms(zone)%surfaceSections(surfSecCounter)%elemNormal(elem,:) = normal
              if(gridDoms(zone)%surfaceSections(surfSecCounter)%isSymmBC) then
                 symDir = gridDoms(zone)%surfaceSections(surfSecCounter)%elemNormal(elem,:)
              end if
              deallocate(points,STAT=ierr)
           end do
           deallocate(center,area,normal,STAT=ierr)
           surfSecCounter = surfSecCounter + 1

        end if
     end do
  end do
!stop
end subroutine getSurfaceElementCenterAndArea
