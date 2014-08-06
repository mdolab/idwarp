subroutine getNPatches(nPatches)

  ! Return the number of wall patches we have on this proc
  use precision
  use constants
  use gridData
 
  implicit none

  ! Output Variables
  integer(kind=intType),intent(out) :: nPatches

  ! Working Variables
  integer(kind=intType) :: zone, sec
  integer(kind=intType) :: surfSecCounter,nSections

  nPatches = 0_intType

  ! loop over the zones and read
  do zone = 1,nZones
     surfSecCounter = 1
     nSections = gridDoms(zone)%nSections 
     do sec=1,nSections
        if (.not. gridDoms(zone)%isVolumeSection(sec)) then
           if(gridDoms(zone)%surfaceSections(surfSecCounter)%isWallBC)then
              nPatches = nPatches + 1
           end if
           surfSecCounter = surfSecCounter + 1
        end if
     end do
  end do
           
end subroutine getNPatches

subroutine getPatchName(iPatch, patchName) 
  
  ! Get names one at a time since f2py can't do arrays of strings
  use precision
  use constants
  use gridData

  implicit none

  ! Input Variables
  integer(kind=intType),intent(in) :: iPatch

  ! Output Variables 
  character*32, intent(out) :: patchName

  ! Working
  integer(kind=intType) :: patchCount
  integer(kind=intType) :: zone,sec
  integer(kind=intType) :: surfSecCounter, nSections
  !print *,'on Entry',iPatch!,patchName
  patchCount = 0
  do zone = 1,nZones
     surfSecCounter = 1
     nSections = gridDoms(zone)%nSections 
     do sec=1,nSections
        if (.not. gridDoms(zone)%isVolumeSection(sec)) then
           if(gridDoms(zone)%surfaceSections(surfSecCounter)%isWallBC)then
              !print *,'matching',patchCount,iPatch
              if (patchCount == iPatch) then
                 patchName = &
                      gridDoms(zone)%surfaceSections(surfSecCounter)%BCFamily
              end if
              patchCount = patchCount + 1
           end if
           surfSecCounter = surfSecCounter + 1
        end if
     end do
  end do

end subroutine getPatchName

subroutine getPatchIndex(iPatch, patchIndex) 
  
  ! Get names one at a time since f2py can't do arrays of strings
  use precision
  use constants
  use gridData

  implicit none

  ! Input Variables
  integer(kind=intType),intent(in) :: iPatch

  ! Output Variables 
  integer(kind=intType),intent(out) ::patchIndex

  ! Working
  integer(kind=intType) :: patchCount
  integer(kind=intType) :: zone,sec
  integer(kind=intType) :: surfSecCounter, nSections
  !print *,'on Entry',iPatch!,patchName
  patchIndex = 0
  patchCount = 0
  do zone = 1,nZones
     surfSecCounter = 1
     nSections = gridDoms(zone)%nSections 
     do sec=1,nSections
        !print *,'index Sections',sec
        if (.not. gridDoms(zone)%isVolumeSection(sec)) then
           if(gridDoms(zone)%surfaceSections(surfSecCounter)%isWallBC)then
              !print *,'matching',patchCount,iPatch,patchCount == iPatch
              if (patchCount == iPatch) then                 
                 patchIndex = surfSecCounter
              end if
              !print *,'patchIndex',patchIndex,surfSecCounter
              patchCount = patchCount + 1
           end if
           surfSecCounter = surfSecCounter + 1
        end if
     end do
  end do
  !print *,'patchIndex',patchIndex,iPatch
  !stop
end subroutine getPatchIndex

! subroutine getPatchSize(iPatch, patchSize)

!   ! Get the patch sizes
!   use precision
!   use constants
!   use gridData

!   implicit none

!   ! Input Variables
!   integer(kind=intType),intent(in) :: iPatch

!   ! Output Variables 
!   integer(kind=intType), intent(out) :: patchSize

!   ! Working
!   integer(kind=intType) :: patchCount
!   integer(kind=intType) :: zone, sec
!   integer(kind=intType) :: surfSecCounter, nSections
!   print *,'ipatch Size',iPatch

!   patchCount = 0
!   do zone = 1,nZones
!      surfSecCounter = 1
!      nSections = gridDoms(zone)%nSections 
!      do sec=1,nSections
!         if (.not. gridDoms(zone)%isVolumeSection(sec)) then
!            if(gridDoms(zone)%surfaceSections(surfSecCounter)%isWallBC)then
!               if (patchCount == iPatch) then
!                  patchSize = &
!                       gridDoms(zone)%surfaceSections(surfSecCounter)%nSurf
!               end if
!               patchCount = patchCount + 1
!            end if
!            surfSecCounter = surfSecCounter + 1
!         end if
!      end do
!   end do
!   print *,'patchsize',patchSize

! end subroutine getPatchSize

subroutine getPatchSize(patchIndex, patchSize)

  ! Get the patch sizes
  use precision
  use constants
  use gridData

  implicit none

  ! Input Variables
  integer(kind=intType),intent(in) :: patchIndex

  ! Output Variables 
  integer(kind=intType), intent(out) :: patchSize

  ! Working
  integer(kind=intType) :: pt, pointPatch
!  print *,'ipatch Size',patchIndex

  ! loop over the unique surface nodes to determine the number of unique nodes
  ! from this surface patch
  patchSize = 0
  do pt = 1,nUniqueSurfPoints
     !base this off of the first connected element
     pointPatch = uniqueSurfaceNodes(pt)%connectedElements(1,2)
     !print *,'pointPatch',pt,pointPatch,patchIndex,'ce', uniqueSurfaceNodes(pt)%connectedElements(:,2)
     if(pointPatch .eq. patchIndex)then
        ! this point is on this patch, append add to counter
        patchSize = patchSize + 1
        !print *,'pointPatch',pt,pointPatch,patchIndex,patchSize
     end if
  end do
  ! print *,'patchsize',patchSize
  ! stop
end subroutine getPatchSize

subroutine getPatchIndexList(patchIndex, patchSize, patchIndices)

  ! Get the patch sizes
  use precision
  use constants
  use gridData

  implicit none

  ! Input Variables
  integer(kind=intType),intent(in) :: patchIndex,patchSize

  ! Output Variables 
  integer(kind=intType),dimension(patchSize), intent(out) :: patchIndices

  ! Working
  integer(kind=intType) :: pt, pointPatch,patchCounter
  !print *,'ipatch Size',patchIndex

  ! loop over the unique surface nodes to determine the number of unique nodes
  ! from this surface patch
  patchCounter = 1
  do pt = 1,nUniqueSurfPoints
     !base this off of the first connected element
     pointPatch = uniqueSurfaceNodes(pt)%connectedElements(1,2)
     if(pointPatch .eq. patchIndex)then
        ! this point is on this patch, append index to list
        patchIndices(patchCounter)=pt
        patchCounter = patchCounter + 1
     end if
  end do

  
end subroutine getPatchIndexList


! subroutine getPatches(points, ndof)

!   ! get the surface points
!   use precision
!   use constants
!   use gridData

!   implicit none
 
!   ! Input/Ouput
!   integer(kind=intType), intent(in) :: ndof
!   real(kind=realType), intent(inout) :: points(ndof)

!   ! Working
!   integer(kind=intType) :: ierr, zone, sec, nSections, ii
!   integer(kind=intType) :: idxStart, idxEnd, glbIdx, pt, surfSecCounter

!   ii = 0
!   do zone = 1,nZones
!      surfSecCounter = 1
!      nSections = gridDoms(zone)%nSections 
!      do sec=1,nSections
!         if (.not. gridDoms(zone)%isVolumeSection(sec)) then
!            if(gridDoms(zone)%surfaceSections(surfSecCounter)%isWallBC)then
!               do pt = 1,gridDoms(zone)%surfaceSections(surfSecCounter)%nSurf
!                  idxStart = ii*3+1
!                  idxEnd = (ii+1)*3
!                  glbIdx = gridDoms(zone)%surfaceSections(surfSecCounter)%uniqueSurfaceNodes(pt)%globalIndex
!                  points(idxStart:idxEnd) = gridDoms(zone)%points(glbIdx,:)
!                  ii = ii + 1
!               end do
!            end if
!            surfSecCounter = surfSecCounter + 1
!         end if
!      end do
!   end do
         
 
! end subroutine getPatches

subroutine getPatches(points, ndof)

  ! get the surface points
  use precision
  use constants
  use gridData

  implicit none
 
  ! Input/Ouput
  integer(kind=intType), intent(in) :: ndof
  real(kind=realType), intent(inout) :: points(ndof)

  ! Working
  integer(kind=intType) :: pt,i
  integer(kind=intType) :: idxStart, idxEnd
  !print *,'check that the number of nodes is correct',ndof,nUniqueSurfPoints,3*nUniqueSurfPoints
  do pt = 1,nUniqueSurfPoints
     do i = 1,3
        idxStart = (pt-1)*3+1
        idxEnd = (pt)*3
        points(idxStart:idxEnd) = uniqueSurfaceNodes(pt)%loc
     end do
  end do
         
 
end subroutine getPatches


subroutine getSurfaceCoordinates(nPoints,indices,points)

  ! get the surface points
  use precision
  use constants
  use gridData

  implicit none
 
  ! Input/Ouput
  integer(kind=intType), intent(in) :: nPoints
  integer(kind=intType),dimension(nPoints), intent(in)::indices
  real(kind=realType), intent(inout) :: points(nPoints*3)

  ! Working
  integer(kind=intType) :: pt,i
  integer(kind=intType) :: idxStart, idxEnd,idx
  !print *,'check that the number of nodes is correct',ndof,nUniqueSurfPoints,3*nUniqueSurfPoints
  do pt = 1,nPoints
     idx = indices(pt)
     do i = 1,3
        idxStart = (pt-1)*3+1
        idxEnd = (pt)*3
        points(idxStart:idxEnd) = uniqueSurfaceNodes(idx)%loc
     end do
  end do
         
end subroutine getSurfaceCoordinates

subroutine setSurfaceCoordinates(nPoints,indices,points)

  ! get the surface points
  use precision
  use constants
  use gridData

  implicit none
 
  ! Input/Ouput
  integer(kind=intType), intent(in) :: nPoints
  integer(kind=intType),dimension(nPoints), intent(in) :: indices
  real(kind=realType),dimension(nPoints*3), intent(in) :: points

  ! Working
  integer(kind=intType) :: pt,i
  integer(kind=intType) :: idxStart, idxEnd,idx
  !print *,'check that the number of nodes is correct',ndof,nUniqueSurfPoints,3*nUniqueSurfPoints
  do pt = 1,nPoints
     idx = indices(pt)
     do i = 1,3
        idxStart = (pt-1)*3+1
        idxEnd = (pt)*3
        uniqueSurfaceNodes(idx)%loc=points(idxStart:idxEnd)
        
     end do
     !print *,'connected elements',uniqueSurfaceNodes(idx)%connectedElements
  end do
         
end subroutine setSurfaceCoordinates


