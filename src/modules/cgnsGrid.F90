module cgnsGrid
  ! 
  ! A module to hold the data structures to store information related
  ! to a structured or unstructured CGNS grid
  !
  use constants
  implicit none
  save

  type bctype
     character(maxCGNSNameLen) :: name
     integer(kind=intType) :: type
     integer(kind=intType) :: ptRange(3, 2)
     character(maxCGNSNameLen) :: family

     ! Required for unstructured grids
     integer(kind=intTYpe) :: nBCElem
     integer(kind=intType) :: nBCNodes
     integer(kind=intType), dimension(:), allocatable :: BCElements
     integer(kind=intType), dimension(:), pointer :: elemPtr, elemConn
     real(kind=realType), dimension(:,:), pointer :: elemNodes
  end type bctype

  type B2B
     ! Only for structured grids
     character(maxCGNSNameLen) :: name
     character(maxCGNSNameLen) :: donorName
     integer(kind=intType) :: ptRange(3, 2)
     integer(kind=intType) :: donorRange(3, 2)
     integer(kind=intType) :: transform(3)
  end type B2B

  type sectionDataType
     ! Only for unstructured grid
     character(maxCGNSNameLen) :: name
     logical :: isSurface
     integer(kind=intType):: nElem
     integer(kind=intType):: elemStart, elemEnd
     integer(kind=intType), dimension(:), pointer :: elemPtr, elemConn
  end type sectionDataType

  type zoneDataType
     integer(kind=intType) :: il, jl, kl
     integer(kind=intType):: nVertices, nElements
     character(maxCGNSNameLen) :: name
     type(bcType), dimension(:), allocatable :: bocos
     type(B2B), dimension(:), allocatable :: B2Bs
     type(sectionDataType), dimension(:), pointer :: sections

  end type zoneDataType

  ! List of the zones
  type(zoneDataType), dimension(:), allocatable :: zones

  ! Deduced information of the wall surfaces
  real(kind=realType), dimension(:), allocatable :: wallPoints
  integer(kind=intType), dimension(:), allocatable :: wallConn
  integer(kind=intType), dimension(:, :), allocatable :: wallSizes
  integer(kind=intType), dimension(:), allocatable :: wallPtr
  integer(kind=intType), dimension(:), allocatable :: wallPatchPtr
  character(maxCGNSNameLen), dimension(:), allocatable :: wallNames

end module cgnsGrid
