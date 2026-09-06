# PETSc Get/Restore audit

A comprehensive sweep of the hand-written IDWarp source for PETSc "borrow"
(`*Get*` / `*Restore*`) calls, checking that every borrow is released and
released against the *correct* object.

## Which PETSc calls need a `Restore`

PETSc splits its accessors into two groups. Only the *borrow* accessors hand
back an internal pointer that the caller must return with a matching
`*Restore*` call; the *copy/scalar* accessors return owned data and have no
restore.

Confirmed against the installed headers (`petsc-3.23.7/include/petscvec.h`):

| Call | Restore counterpart? |
|------|----------------------|
| `VecGetArray` | **yes** → `VecRestoreArray` |
| `VecGetArrayRead`, `VecGetArrayF90`, `ISGetIndices`, `MatDenseGetArray`, `MatGetRow`, `VecGhostGetLocalForm`, `VecGetSubVector`, … | yes (not used in IDWarp) |
| `VecGetOwnershipRange` | no (returns scalars) |
| `VecGetOwnershipRanges` | no (returns a pointer *into* the vector's `PetscLayout`; must be copied, never restored/freed) |
| `VecGetSize` | no (scalar) |
| `VecGetValues` | no (copies values into a caller array) |

The only borrow accessor actually used in the hand-written source is
`VecGetArray`. The generated adjoint code (`src/adjoint/outputForward`,
`src/adjoint/outputReverse`) and the f2py wrapper (`src/f2py`) contain no
borrow accessors, so nothing there needs a restore.

## Result

20 `VecGetArray` calls, 20 `VecRestoreArray` calls. Counts match, and after
the earlier fixes on this branch every borrow is now released — but **two
releases target the wrong object**. A `VecRestoreArray` against the wrong
Vec/pointer means the array that *was* borrowed is never really returned,
while an object that was never borrowed is "restored" instead.

### Finding 1 — `getCommonVolumeCoordinates.F90`: restores the wrong Vec

```fortran
call VecGetArray(commonGridVec, xx, ierr)   ! line 23: borrow from commonGridVec
...
call VecRestoreArray(Xv, xx, ierr)          ! line 29: restores Xv (WRONG)
```

The array is borrowed from `commonGridVec` but restored against `Xv`.
(The sibling routine `setCommonVolumeCoordinates` in the same file restores
`commonGridVec` correctly.)

- **Impact — correctness/robustness.** `commonGridVec`'s borrow is never
  released, and `Xv` — which was never borrowed here — has its object state
  bumped spuriously (invalidating any cached norms). In a debug PETSc build
  the mismatched pointer/Vec can trip an internal check and abort through
  `EChk`. This is a read path (`gridNodes = xx`), so no data is corrupted,
  but the call is unambiguously wrong.
- **Fix:** restore `commonGridVec`.

### Finding 2 — `warpDerivFwd.F90`: restores the wrong pointer

```fortran
call VecGetArray(dXsLocal, XsPtrd, ierr)    ! line 51: borrow into XsPtrd
...
call VecRestoreArray(dXsLocal, dXsPtr, ierr) ! line 105: restores dXsPtr (WRONG)
```

`XsPtrd` and `dXsPtr` are two distinct module pointers declared in
`gridData.F90`. The array is borrowed into `XsPtrd` but restored via
`dXsPtr`, which was never borrowed in this routine.

- **Impact — potential crash / undefined behaviour.** `XsPtrd`'s borrow is
  never released, and `dXsPtr` is passed to `VecRestoreArray` while
  unassociated or stale, which is undefined (PETSc error or segfault).
  Guarded by `#ifndef USE_COMPLEX`, so only the real build is affected.
- **Fix:** restore `XsPtrd`.

## Adjacent observation (not a restore issue, not fixed)

`warpDerivFwd.F90:61` reads `nVol = size(XvPtr) / 3`, but `XvPtr` is never
borrowed in that routine (it uses `Xv0Ptr`). This is a stale-module-pointer
read, out of scope for this restore-pairing audit — flagged for a follow-up.

## Fixes applied

Both wrong-target restores corrected to match their borrow.

## Verification

Clean builds and the full `testflo` suite were run for real **and** complex,
against PETSc 3.23.7 and 3.25.5 (both `real-opt` and `complex-opt` arches):

- All four builds succeed.
- 12 tests pass; 8 fail. The 8 failures are exclusively the **parallel
  (`N_PROCS=2`) `test_USMesh` variants** (`test_comesh`, `test_inflate_cube`,
  `test_omesh`, `test_sym_mesh`, real and complex). They are a numeric
  mismatch (~1.3e-5 rel. on "Sum of dxs"), identical across both PETSc
  versions.
- These failures **pre-date this work**: rebuilding with the restore fixes
  stashed reproduces the same `test_sym_mesh` failure. The restore fixes are
  numerically inert (data is copied before the restore; `VecRestoreArray`
  only nulls the pointer and bumps object state), so they neither cause nor
  cure the parallel mismatch. All serial `test_USMesh` cases and all
  `test_MultiUSMesh` cases (serial and parallel) pass.
