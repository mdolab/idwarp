! ====================================================================
! File: createSections.f90
! Author: C.A.(Sandy) Mader
! Date Started: July 30, 2014
! Date Modified:

subroutine createSections(zone,nSections)
  use precision
  use gridData
  implicit none

  ! Subroutine variables
  integer(kind=intType),intent(in)::zone,nSections

  ! Local variables
  integer(kind=intType)::ierr
  
  ! begin execution
  gridDoms(zone)%nSections = nSections

  ! allocate the logical to track the section type
  allocate(gridDoms(zone)%isVolumeSection(nSections),STAT=ierr)
  gridDoms(zone)%isVolumeSection=.True.

  ! initialize the family counters and list
  nwallFamilies = 0
  familyList(:) = ''


end subroutine createSections


subroutine setSectionTypes(zone,nSection,isVolumeSec)
  use precision
  use gridData
  implicit none

  ! Subroutine variables
  integer(kind=intType),intent(in)::zone,nSection
  logical, dimension(nSection),intent(in):: isVolumeSec
  ! Local variables
  integer(kind=intType)::ierr,nVolSec,nSurfSec,sec
  
  ! begin execution
  nVolSec = 0
  nSurfSec = 0
  do sec=1,nSection
     if(isVolumeSec(sec))then
        nVolSec = nVolSec+1
     else
        nSurfSec = nSurfSec+1
     end if
  end do
  gridDoms(zone)%nVolSec = nVolSec
  gridDoms(zone)%nSurfSec = nSurfSec
  allocate(gridDoms(zone)%volumeSections(gridDoms(zone)%nVolSec),STAT=ierr)
  allocate(gridDoms(zone)%surfaceSections(gridDoms(zone)%nSurfSec),STAT=ierr)
  !Fix isVolumeSec(sec)! pass in an array, and count numbers from that.
  gridDoms(zone)%isVolumeSection(:) = isVolumeSec(:)
end subroutine setSectionTypes

subroutine setSectionData(zone,sec,secCounter,nElem,nNodeIdx,secName,eBeg,eEnd,&
     nodeIndexList,elemPtr)
  use precision
  use gridData
  implicit none

  ! Subroutine variables
  integer(kind=intType),intent(in)::zone,sec,secCounter,nElem,nNodeIdx
  integer(kind=intType),intent(in)::eBeg,eEnd
  character*32 :: secName
  integer(kind=intType),dimension(nNodeIdx)::nodeIndexList
  integer(kind=intType),dimension(nElem+1):: elemPtr

  ! local variables
  integer(kind=intType):: ierr,surfSecCounter
  ! begin execution
  if (gridDoms(zone)%isVolumeSection(sec)) then
     ! this is a group of volume elements
     allocate(gridDoms(zone)%volumeSections(secCounter)%elemPtr(nElem+1),STAT=ierr)
     allocate(gridDoms(zone)%volumeSections(secCounter)%elemConn(nNodeIdx),STAT=ierr)
     gridDoms(zone)%volumeSections(secCounter)%elemPtr = elemPtr
     gridDoms(zone)%volumeSections(secCounter)%elemConn = nodeIndexList
     gridDoms(zone)%volumeSections(secCounter)%nElem = nElem
     gridDoms(zone)%volumeSections(secCounter)%secName = secName
     !gridDoms(zone)%volumeSections(secCounter)%elemType = type
     !gridDoms(zone)%volumeSections(secCounter)%nConn = nConn 
     gridDoms(zone)%volumeSections(secCounter)%elemStart = eBeg
     gridDoms(zone)%volumeSections(secCounter)%elemEnd = eEnd
  else
     ! this is a group of surface elements
     allocate(gridDoms(zone)%surfaceSections(secCounter)%elemPtr(nElem+1),STAT=ierr)
     allocate(gridDoms(zone)%surfaceSections(secCounter)%elemConn(nNodeIdx),STAT=ierr)

     allocate(gridDoms(zone)%surfaceSections(secCounter)%elemCenter(nElem,3),STAT=ierr)
     allocate(gridDoms(zone)%surfaceSections(secCounter)%elemArea(nElem,3),STAT=ierr)
     allocate(gridDoms(zone)%surfaceSections(secCounter)%elemNormal(nElem,3),STAT=ierr)
     allocate(gridDoms(zone)%surfaceSections(secCounter)%elemAreaMag(nElem),STAT=ierr)
     !print *,'allocated',shape(gridDoms(zone)%surfaceSections(secCounter)%elemNormal),nElem
     gridDoms(zone)%surfaceSections(secCounter)%elemPtr = elemPtr
     gridDoms(zone)%surfaceSections(secCounter)%elemConn = nodeIndexList
     gridDoms(zone)%surfaceSections(secCounter)%nElem = nElem
     gridDoms(zone)%surfaceSections(secCounter)%secName = secName
     !gridDoms(zone)%surfaceSections(secCounter)%elemType = type
     !gridDoms(zone)%surfaceSections(secCounter)%nConn = nConn 
     gridDoms(zone)%surfaceSections(secCounter)%elemStart = eBeg
     gridDoms(zone)%surfaceSections(secCounter)%elemEnd = eEnd        
  end if

end subroutine setSectionData

subroutine setBoundaryData(zone,sec,surfSecCounter,famName,isWall,isSymm,&
     isBoundary)
  use precision
  use gridData
  implicit none

  ! Subroutine Variables
  integer(kind=intType),intent(in)::zone,sec,surfSecCounter
  character*32::famName
  logical:: isWall,isSymm,isBoundary

  ! Local Variables
  integer(kind=intType):: famID

 ! We can also generate the wall family list
  
  if (.not. gridDoms(zone)%isVolumeSection(sec)) then
     gridDoms(zone)%surfaceSections(surfSecCounter)%isWallBC=.False.
     gridDoms(zone)%surfaceSections(surfSecCounter)%isBoundaryBC=.False.
     gridDoms(zone)%surfaceSections(surfSecCounter)%isSymmBC=.False.
     !print *,'boundaries',sec,surfSecCounter,isWall,isSymm,isBoundary,famName, gridDoms(zone)%surfaceSections(surfSecCounter)%secName, gridDoms(zone)%surfaceSections(surfSecCounter)%elemStart
     if (isWall) then
        !print *,'is Wall',zone,surfSecCounter
        gridDoms(zone)%surfaceSections(surfSecCounter)%isWallBC=.True.
        gridDoms(zone)%surfaceSections(surfSecCounter)%BCFamily= famName
        !this is a wall family, add to the family list
        call checkInFamilyList(familyList, famName, nwallFamilies, famID)
        !print *,'famID',famID,nwallFamilies
        if (famID == 0) then
           nwallFamilies = nwallFamilies + 1
           familyList(nwallFamilies) = famName
           !print *,'nfamilies',nwallFamilies,familyList(nwallFamilies)
        end if
     elseif(isBoundary)then
        !print *,'is farfield',zone,surfSecCounter
        gridDoms(zone)%surfaceSections(surfSecCounter)%isBoundaryBC=.True.
     elseif(isSymm) then
        !print *,' is symm',zone,surfSecCounter
        gridDoms(zone)%surfaceSections(surfSecCounter)%isSymmBC = .True.
        hasSymmetry = .True.
     else
        print *,'Unrecongnized boundary Type... exiting'
        stop
     end if
  end if

end subroutine setBoundaryData
