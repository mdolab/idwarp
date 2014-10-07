module gridData
  ! 
  ! A module to hold the data structures for the grid data
  !

  use constants
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
  Vec commonGridVec, commonSymmMask
  Vec solverGridVec
  VecScatter common_to_solver
  VecScatter common_to_warp
#endif

  ! Pointers into the grid vecs
  real(kind=realType), pointer, dimension(:) :: XsPtr, dXsPtr
  real(kind=realType), pointer, dimension(:) :: Xv0Ptr, Xvptr
  integer(kind=intType) :: warpMeshDOF
  integer(kind=intType) :: commonMeshDOF

#ifndef USE_TAPENADE
  real(kind=realType), pointer, dimension(:) :: XsPtrb
  real(kind=realType), pointer, dimension(:) :: XvPtrb
#endif

  logical commonGridVecSet 
  logical gridIndicesSet 

  character*32, dimension(maxFamilies) :: familyList
  integer(kind=intType) :: nwallFamilies

  real(kind=realType), dimension(:), allocatable :: denomenator, denomenator0
  real(kind=realType), dimension(:, :), allocatable :: numerator

end module gridData
