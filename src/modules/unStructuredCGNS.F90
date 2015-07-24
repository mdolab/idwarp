module unStructuredCGNSGrid
  ! 
  ! A module to hold the data structures to store information related
  ! to a structured CGNS grid
  !
  use constants
  implicit none
  save
 
  type bcType
     integer(kind=intType),dimension(:),allocatable::BCElements
     integer(kind=intType)::nBCElem,BCType
     character*32 :: BCName, famName
  end type bcType

  type sectionDataType
     character*32 ::secName
     integer(kind=intType):: nElem,nConn
     integer(kind=intType):: elemStart, elemEnd
     integer(kind=intType):: elemType
     real(kind=realType), dimension(:), allocatable :: elemPtr,elemConn
     real(kind=realType), dimension(:,:), allocatable :: elemCenter
     real(kind=realType), dimension(:,:), allocatable :: elemArea
     real(kind=realType), dimension(:,:), allocatable :: elemNormal
     real(kind=realType), dimension(:), allocatable :: elemAreaMag
     integer(kind=intType)::nSurf

     ! Section BC Information
     character*32 :: BCFamily
     integer(kind=intType) :: BCType
     logical :: isWallBC,isBoundaryBC,isSymmBC
  end type sectionDataType

  type zoneDataType
     !
     ! nSections:  Number of element sections in this zone. Sections can be 
     !             sets of either surface or volume elements
     ! nVolSec:    Number of volume element sections in this zone
     ! nSurfSec:   Number of surface element sections in this zone
     ! points      Full array of points for this zone. this will be used for
     !             all of the other element definitions
     ! volumeSections: This is a data structure to store the volume related
     !                 section data.
     ! surfaceSections: This is a data structure to store the surface related
     !                  section data.
     ! isVolumeSection: An array of logicals to help sort the section data


     character*32 :: zoneName
     integer(kind=intType)::nSections, nVertices,nElements
     integer(kind=intType)::nVolSec,nSurfSec
     integer(kind=intType)::volSecCounter, surfSecCounter

     type(sectionDataType), dimension(:),allocatable :: volumeSections
     type(sectionDataType), dimension(:),allocatable :: surfaceSections
     logical, dimension(:), allocatable :: isVolumeSection
     type(bcType), dimension(:),allocatable :: BCInfo
     integer(kind=intType) :: nBocos

  end type zoneDataType

  integer(kind=intType):: nZones
  type(zoneDataType), dimension(:), allocatable :: zones

end module unStructuredCGNSGrid
