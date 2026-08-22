module gridData
    !
    ! A module to hold the data structures for the grid data
    !

    use constants

#ifndef USE_TAPENADE

#include <petsc/finclude/petsc.h>

    ! On PETSc 3.22 and older, the post-3.23 Fortran call forms bind silently to the old F77 stubs
    ! and corrupt the vector layout at runtime instead of failing to link, so this floor must be
    ! enforced here at compile time, not discovered later by a mysterious runtime failure.
#if PETSC_VERSION_LT(3,23,0)
#error "IDWarp requires PETSc 3.23 or newer"
#endif

    use petsc
    implicit none

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
    real(kind=realType), dimension(:), allocatable :: denominator, denominator0
    ! target is required: numerator is rank-remapped to the pointer numerator1D for VecPlaceArray,
    ! and a pointer target must be declared target.
    real(kind=realType), dimension(:, :), allocatable, target :: numerator

    ! Symmetry Information
    integer(kind=intType) :: nLoop
    real(kind=realType), dimension(:, :), allocatable :: symmPts, symmNormals

end module gridData

module plot3dSurface

    use constants
    implicit none
    save

    real(kind=realType), dimension(:, :), allocatable :: pts
    integer(kind=intType), dimension(:, :), allocatable :: conn
end module plot3dSurface
