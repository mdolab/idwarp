module gridData
  ! 
  ! A module to hold the data structures for the grid data read in
  ! from files for either CGNS Data or OpenFOAM data
  !

  use precision
  use constants
  use gridTypes
  implicit none
  save

  integer(kind=intType):: nZones,physDim,cellDim
  type(zoneDataType) ,dimension(:), allocatable ::gridDoms
  real(kind=realType):: aExp,bExp,Ldef,alpha

  ! Family Data
  character*32, dimension(maxFamilies) :: familyList
  integer(kind=intType) :: nwallFamilies

  ! All of the wall nodes
  type(surfacePointType),dimension(:),allocatable :: uniqueSurfaceNodes
  integer(kind=intType)::nSurfNodes,nUniqueSurfPoints
  ! Boundary nodes are all of the non-surface, non-symmetry boundary nodes
  type(surfacePointType),dimension(:),allocatable :: uniqueBoundaryNodes
  integer(kind=intType)::nBoundaryNodes,nUniqueBoundaryPoints
  logical:: hasSymmetry
  real(kind=realType),dimension(3)::symDir !Symmetry direction
 
end module gridData
