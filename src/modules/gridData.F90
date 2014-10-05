module gridData
  ! 
  ! A module to hold the data structures for the grid data
  !

  use constants
  use kd_tree
  implicit none
  save
#ifndef USE_TAPENADE
#include "finclude/petsc.h"
#include "finclude/petscvec.h90"

  ! Volume grid vecs
  Vec Xv, Xv0, dXv, XvLocal

  ! Surface Grid vecs
  Vec Xs, dXs, XsLocal

  ! Scatter context going from partitioned Xs to (full) local Copy
  VecScatter XsToXsLocal

  ! Temporary scatter from Xv to full local version
  VecScatter Xvtolocal

  ! Generic index sets
  IS IS1, IS2

  ! Scatter/ (empty) vectors for doing external solver communication
  Vec commonGridVec
  Vec solverGridVec
  VecScatter common_to_solver
  VecScatter common_to_warp
#endif

  ! Pointers into the grid vects
  real(kind=realType), pointer, dimension(:) :: XsPtr
  real(kind=realType), pointer, dimension(:) :: dXsPtr
  real(kind=realType), pointer, dimension(:) :: Xv0Ptr, Xvptr
  integer(kind=intType) :: warpMeshDOF
  integer(kind=intType) :: commonMeshDOF

#ifndef USE_TAPENADE
  real(kind=realType), pointer, dimension(:) :: XsPtrb
  real(kind=realType), pointer, dimension(:) :: XvPtrb
#endif

  logical commonGridVecSet 
  logical gridIndicesSet 
  integer(kind=intType), allocatable, dimension(:) :: facePtr, faceConn
  integer(kind=intType) nFace, lenFaceConn


  character*32, dimension(maxFamilies) :: familyList
  integer(kind=intType) :: nwallFamilies

  integer(kind=intType) :: nUnique
  real(kind=realType) :: Ldef0
  real(kind=realType), dimension(:, :), pointer :: Xu, Xu0
  integer(kind=intType), dimension(:), allocatable :: XuInd
  real(kind=realType), dimension(:, :), allocatable :: XuFact
  real(kind=realType), dimension(:, :, :), pointer :: Mi
  real(kind=realType), dimension(:, :), pointer :: Bi
  real(kind=realType), dimension(:, :), allocatable :: normals0, normals
  real(kind=realType), dimension(:), pointer :: Ai
  real(kind=realType), dimension(:), allocatable :: denomenator, denomenator0
  real(kind=realType), dimension(:, :), allocatable :: numerator

#ifndef USE_TAPENADE
  real(kind=realType), dimension(:, :), allocatable :: Xub
  real(kind=realType), dimension(:, :, :), allocatable :: Mib
  real(kind=realType), dimension(:, :), allocatable :: Bib
  real(kind=realType), dimension(:, :), allocatable :: normals0b, normalsb
  real(kind=realType), dimension(:), allocatable :: Aib
  real(kind=realType), dimension(:), allocatable :: denomenatorb
  real(kind=realType), dimension(:, :), allocatable :: numeratorb

#endif  
  type(tree_master_record), pointer :: mytree
  integer(kind=intType), dimension(:, :), allocatable :: nodeToElem
end module gridData
