module structuredCGNSGrid
  ! 
  ! A module to hold the data structures to store information related
  ! to a structured CGNS grid
  !
  use constants
  implicit none
  save

  type bctype
     character(maxCGNSNameLen) :: name
     integer(kind=intType) :: bocoType
     integer(kind=intType) :: ptRange(3, 2)
     character(maxCGNSNameLen) :: family
  end type bctype

  type B2B
     character(maxCGNSNameLen) :: name
     character(maxCGNSNameLen) :: donorName
     integer(kind=intType) :: ptRange(3, 2)
     integer(kind=intType) :: donorRange(3, 2)
     integer(kind=intType) :: transform(3)
  end type B2B

  type blocktype
     integer :: il, jl, kl
     character(maxCGNSNameLen) :: name
     type(bcType), dimension(:), allocatable :: bocos
     type(B2B), dimension(:), allocatable :: B2Bs
  end type blocktype

  type(blockType), dimension(:), allocatable :: blocks

end module structuredCGNSGrid
