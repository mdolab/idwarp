module gridData
  ! 
  ! A module to hold the data structures for the grid data
  !

  use constants
  implicit none
  save

#ifndef USE_TAPENADE

#include "include/petscversion.h"
#if PETSC_VERSION_MINOR > 5
#include "petsc/finclude/petsc.h"
#include "petsc/finclude/petscvec.h90"
#else
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h90"
#endif

  ! Volume grid vecs
  Vec Xv, Xv0, dXv, XvLocal

  ! Surface Grid vecs
  Vec Xs, dXs, XsLocal, dXsLocal

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
  VecScatter common_to_dxs
#endif

  ! Pointers into the grid vecs
  real(kind=realType), pointer, dimension(:) :: XsPtr, dXsPtr
  real(kind=realType), pointer, dimension(:) :: Xv0Ptr, Xvptr

#ifndef USE_TAPENADE
  real(kind=realType), pointer, dimension(:) :: XsPtrb, XsPtrd
  real(kind=realType), pointer, dimension(:) :: XvPtrb, XVPtrd
#endif

  ! Sizes of the three different mesh sizes:
  integer(kind=intType) :: warpMeshDOF
  integer(kind=intType) :: commonMeshDOF
  integer(kind=intType) :: solverMeshDOF

  ! Logicals determine what is allocated:
  integer(kind=intTYpe) :: gridIndicesSet = 0
  integer(kind=intType) :: commonGridVecSet = 0
  integer(kind=intTYpe) :: initializationSet = 0

  real(kind=realType), dimension(:), allocatable :: d2wall
  real(kind=realType), dimension(:), allocatable :: denomenator, denomenator0
  real(kind=realType), dimension(:, :), allocatable :: numerator
  integer(kind=intType), dimension(:), allocatable :: surfaceIndices

  ! Symmetry Information
  integer(kind=intType) :: nLoop
  real(kind=realType), dimension(:, :), allocatable :: symmPts, symmNormals

end module gridData
