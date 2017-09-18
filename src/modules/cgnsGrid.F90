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
     type(sectionDataType), dimension(:), pointer :: sections

  end type zoneDataType

  ! List of the zones
  type(zoneDataType), dimension(:), allocatable :: zones

  ! Deduced information of the surfaces
  real(kind=realType), dimension(:), allocatable :: surfacePoints
  logical, dimension(:), allocatable :: surfaceIsWall
  logical, dimension(:), allocatable :: surfaceIsSymm
  integer(kind=intType), dimension(:), allocatable :: surfaceConn
  integer(kind=intType), dimension(:, :), allocatable :: surfaceSizes
  integer(kind=intType), dimension(:), allocatable :: surfacePtr
  integer(kind=intType), dimension(:), allocatable :: surfacePatchPtr
  character(maxCGNSNameLen), dimension(:), allocatable :: surfaceNames

  ! Flag if it is structured or not
  logical :: cgnsStructured

  ! List of default names for surfaces if not provided
  character(maxCGNSNameLen), dimension(25) :: defaultFamName

contains

  subroutine getnpatch(npatch)

    implicit none

    integer(kind=inttype), intent(out) :: npatch
    if (allocated(surfaceNames)) then
       npatch = size(surfaceNames)
    else
       npatch = 0
    end if

  end subroutine getnpatch

  subroutine getsurf(i, surf)

    implicit none

    character(len=maxCGNSNameLen), intent(out) :: surf
    integer(kind=intType), intent(in) :: i

    if (allocated(surfaceNames)) then
       surf = surfaceNames(i)
    end if

  end subroutine getsurf

end module cgnsGrid
