module gridData
  ! 
  ! A module to hold the data structures for the grid data
  !

  use constants
  implicit none
  save
#include "finclude/petsc.h"
#include "finclude/petscvec.h90"

  ! Volume grid vecs
  Vec Xv, Xv0, dXv, XvLocal
  integer(kind=intType) :: warpMeshDOF

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

  logical commonGridVecSet 
  logical gridIndicesSet 
  integer(kind=intType), allocatable, dimension(:) :: facePtr, faceConn
  integer(kind=intType) nFace, lenFaceConn

  character*32, dimension(maxFamilies) :: familyList
  integer(kind=intType) :: nwallFamilies

  integer(kind=intType) :: nUnique
  real(kind=realType) :: Ldef0
  real(kind=realType), dimension(:, :), allocatable :: Xu, Xu0
  integer(kind=intType), dimension(:), allocatable :: XuInd
  real(kind=realType), dimension(:, :), allocatable :: XuFact
  real(kind=realType), dimension(:, :, :), allocatable :: Mi
  real(kind=realType), dimension(:, :), allocatable :: Bi
  real(kind=realType), dimension(:, :), allocatable :: normals0, normals
  real(kind=realType), dimension(:), allocatable :: Ai
  integer(kind=intType), dimension(:, :), allocatable :: nodeToElem
end module gridData
