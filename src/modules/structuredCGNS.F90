module structuredCGNSGrid
  ! 
  ! A module to hold the data structures to store information related
  ! to a structured CGNS grid
  !
  use constants
  implicit none
  save

  type bctype
     character*32 name
     integer(kind=intType) :: bocoType
     integer(kind=intType) :: ptRange(3, 2)
     character*32 family
  end type bctype

  type B2B
     character*32 name
     character*32 donorName
     integer(kind=intType) :: ptRange(3, 2)
     integer(kind=intType) :: donorRange(3, 2)
     integer(kind=intType) :: transform(3)
  end type B2B

  type blocktype
     integer :: il, jl, kl
     character*32 :: name
     type(bcType), dimension(:), allocatable :: bocos
     type(B2B), dimension(:), allocatable :: B2Bs
  end type blocktype

  type(blockType), dimension(:), allocatable :: blocks

end module structuredCGNSGrid
