module unStructuredCGNSGrid
  ! 
  ! A module to hold the data structures to store information related
  ! to a structured CGNS grid
  !
  use constants
  implicit none
  save
 
  ! type sectionDataType
     
  !    character*32 ::secName
  !    integer(kind=intType):: nElem,nConn
  !    integer(kind=intType):: elemType
  !    integer(kind=intType), dimension(:,:), allocatable :: elements

  ! end type sectionDataType


  ! type zoneDataType
  !    !
  !    ! nSections:  Number of element sections in this zone. Sections can be 
  !    !             sets of either surface or volume elements
  !    ! nVolSec:    Number of volume element sections in this zone
  !    ! nSurfSec:   Number of surface element sections in this zone
  !    ! points      Full array of points for this zone. this will be used for
  !    !             all of the other element definitions
  !    ! volumeSections: This is a data structure to store the volume related
  !    !                 section data.
  !    ! surfaceSections: This is a data structure to store the surface related
  !    !                  section data.
  !    ! isVolumeSection: An array of logicals to help sort the section data


  !    integer(kind=intType)::nSections
  !    integer(kind=intType)::nVolSec,nSurfSec
  !    integer(kind=intType)::volSecCounter, surfSecCounter

  !    real(kind=realType), dimension(:,:),allocatable :: points
     
  !    type(sectionDataType), dimension(:),allocatable :: volumeSections
  !    type(sectionDataType), dimension(:),allocatable :: surfaceSections
  !    logical,dimension(:),allocatable :: isVolumeSection

  ! end type zoneDataType
 ! type elementDataType
     
  !    integer(kind=intType):: elemType
  !    integer(kind=intType),dimension(:),allocatable::pointList

  ! end type elementDataType
  type bcType
     integer(kind=intType),dimension(:),allocatable::BCElements
     integer(kind=intType)::nBCElem,BCType
     character*32 :: BCName, famName
  end type bcType
  type surfacePointType
     ! Update copySurfacePoint call if this type is changed
     integer(kind=intType):: globalIndex
     ! integer(kind=intType),dimension(:),allocatable::connectedNodes
     ! integer(kind=intType),dimension(:),allocatable::connectedElements
     !integer(kind=intType),dimension(10)::connectedNodes
     integer(kind=intType):: proc, zone!, surfSection
     integer(kind=intType),dimension(30,2)::connectedElements
     real(kind=realType),dimension(3)::loc,loc0
     real(kind=realType),dimension(3)::normal, normal0
     real(kind=realType),dimension(3,3)::Mi
     real(kind=realType),dimension(3)::bi
     real(kind=realType):: Ai,Wi
     real(kind=realType),dimension(3)::Si

     real(kind=realType),dimension(3)::symloc,symloc0
     real(kind=realType),dimension(3)::symNormal, symNormal0
     real(kind=realType),dimension(3,3)::symMi
     real(kind=realType),dimension(3)::symbi
     real(kind=realType):: symAi,symWi
     real(kind=realType),dimension(3)::symSi

  end type surfacePointType

  type sectionDataType
     
     character*32 ::secName
     integer(kind=intType):: nElem,nConn
     integer(kind=intType):: elemStart, elemEnd
     integer(kind=intType):: elemType
     !integer(kind=intType), dimension(:,:), allocatable :: elements
     real(kind=realType), dimension(:), allocatable :: elemPtr,elemConn
     real(kind=realType), dimension(:,:), allocatable :: elemCenter
     real(kind=realType), dimension(:,:), allocatable :: elemArea
     real(kind=realType), dimension(:,:), allocatable :: elemNormal
     real(kind=realType), dimension(:), allocatable :: elemAreaMag

     type(surfacePointType),dimension(:),allocatable :: uniqueSurfaceNodes
     integer(kind=intType)::nSurf

     ! Section BC Information
     character*32 :: BCFamily
     integer(kind=intType):: BCType
     logical :: isWallBC,isBoundaryBC,isSymmBC

     ! nodalAreas
     ! nodalWi
     ! nodaldisplacement
     ! nodalRotation
     ! nodalSi
     

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

     real(kind=realType), dimension(:,:),allocatable :: points,points0
     real(kind=realType), dimension(:,:),allocatable :: dx
     logical, dimension(:),allocatable ::isSurfaceNode
     
     type(sectionDataType), dimension(:),allocatable :: volumeSections
     type(sectionDataType), dimension(:),allocatable :: surfaceSections
     logical,dimension(:),allocatable :: isVolumeSection
     type(bcType),dimension(:),allocatable :: BCInfo
     integer(kind=intType)::nBocos

  end type zoneDataType


  integer(kind=intType):: nZones
  type(zoneDataType) ,dimension(:), allocatable ::zones
  logical:: hasSymmetry !do we still need this?
end module unStructuredCGNSGrid
