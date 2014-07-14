module gridData
  ! 
  ! A module to hold the data structures for the grid data read in
  ! from files for either CGNS Data or OpenFOAM data
  !

  use precision
  implicit none
  save

  ! type elementDataType
     
  !    integer(kind=intType):: elemType
  !    integer(kind=intType),dimension(:),allocatable::pointList

  ! end type elementDataType

  type sectionDataType
     
     character*32 ::secName
     integer(kind=intType):: nElem,nConn
     integer(kind=intType):: elemType
     integer(kind=intType), dimension(:,:), allocatable :: elements

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


     integer(kind=intType)::nSections
     integer(kind=intType)::nVolSec,nSurfSec
     integer(kind=intType)::volSecCounter, surfSecCounter

     real(kind=realType), dimension(:,:),allocatable :: points
     
     type(sectionDataType), dimension(:),allocatable :: volumeSections
     type(sectionDataType), dimension(:),allocatable :: surfaceSections
     logical,dimension(:),allocatable :: isVolumeSection

  end type zoneDataType


  integer(kind=intType):: nZones
  type(zoneDataType) ,dimension(:), allocatable ::gridDoms


end module gridData
