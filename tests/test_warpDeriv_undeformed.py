# =============================================================================
# Standard Python modules
# =============================================================================
import os
import copy

# =============================================================================
# External Python modules
# =============================================================================
import numpy
import unittest
from mpi4py import MPI
from parameterized import parameterized_class

# =============================================================================
# Extension modules
# =============================================================================
from idwarp import USMesh

baseDir = os.path.dirname(os.path.abspath(__file__))  # Path to current folder


def check_warpDeriv_at_undeformed_mesh(testcase, meshOptions, h=1e-6, tol=1e-5, seed=314):
    """Check both warping derivatives against a central finite difference of
    ``warpMesh`` itself, linearized at the UNDEFORMED mesh.

    Why this test exists
    --------------------
    At an undeformed mesh every surface node has current normal == initial
    normal.  ``getRotationMatrix3d`` used to build the nodal rotation from an
    axis-angle form guarded by ``if (axisMag < sqrt(eps))``, and the Tapenade
    derivatives of that branch were identically zero -- so ``warpDeriv`` and
    ``warpDerivFwd`` silently dropped the entire rotation contribution whenever
    they were linearized at a baseline configuration, which is exactly where
    the first optimization iterate and every gradient check are evaluated.

    Why it is written this way
    --------------------------
    The reference is a finite difference of the *primal*, not the other AD mode.
    A forward-vs-reverse dot-product test cannot detect this class of defect:
    both modes were dead on the same branch, so they agreed with each other to
    machine precision while both were wrong.  Deforming the mesh before
    differentiating also hides it, because the old code was correct once the
    normals had rotated away from their initial values.  So: no deformation,
    and the primal is the reference.

    With the defect present this check reads a relative error of O(1) on both
    modes (the AD response norm is roughly 1/8 of the true one on ``o_mesh``).
    With the fix it reads ~1e-7 at h = 1e-6, and the error scales as h^2.
    """
    comm = MPI.COMM_WORLD
    mesh = USMesh(options=meshOptions, comm=comm)

    # --- Linearize at the undeformed configuration: do NOT deform first ---
    Xs0 = mesh.getSurfaceCoordinates().copy()
    mesh.setSurfaceCoordinates(Xs0)
    mesh.warpMesh()

    rng = numpy.random.default_rng(seed)
    dXs = rng.random(Xs0.shape) - 0.5  # generic direction: rotates the nodal normals

    # --- Forward-mode AD ---
    dXv_ad = mesh.warpDerivFwd(dXs, solverVec=False).copy()

    # --- Central finite difference of the primal, same direction ---
    mesh.setSurfaceCoordinates(Xs0 + h * dXs)
    mesh.warpMesh()
    Xv_p = mesh.getWarpGrid().copy()
    mesh.setSurfaceCoordinates(Xs0 - h * dXs)
    mesh.warpMesh()
    Xv_m = mesh.getWarpGrid().copy()
    dXv_fd = (Xv_p - Xv_m) / (2.0 * h)

    # --- Back to the undeformed mesh before the reverse-mode call ---
    mesh.setSurfaceCoordinates(Xs0)
    mesh.warpMesh()

    # --- Reverse-mode AD, checked against the SAME finite difference through
    #     the identity <dXv, dXv/dXs . dXs> == <(dXv/dXs)^T . dXv, dXs> ---
    dXv = rng.random(dXv_fd.shape) - 0.5
    mesh.warpDeriv(dXv, solverVec=False)
    dXs_ad = mesh.getdXs().copy()

    # --- Global norms (the vectors are distributed) ---
    num = comm.allreduce(numpy.sum((dXv_ad - dXv_fd) ** 2), op=MPI.SUM)
    den = comm.allreduce(numpy.sum(dXv_fd**2), op=MPI.SUM)
    fwd_relErr = numpy.sqrt(num / den)

    lhs = comm.allreduce(numpy.sum(dXv * dXv_fd), op=MPI.SUM)
    rhs = comm.allreduce(numpy.sum(dXs_ad * dXs), op=MPI.SUM)
    rev_relErr = abs(lhs - rhs) / abs(lhs)

    testcase.assertLess(
        fwd_relErr,
        tol,
        msg=f"warpDerivFwd disagrees with a central FD of warpMesh at the undeformed mesh: rel err {fwd_relErr:.3e}",
    )
    testcase.assertLess(
        rev_relErr,
        tol,
        msg=f"warpDeriv disagrees with a central FD of warpMesh at the undeformed mesh: rel err {rev_relErr:.3e}",
    )


# The derivative routines exist only in the real build, so there is no complex variant.
test_params = [
    {"N_PROCS": 1},
    {"N_PROCS": 2, "name": "parallel"},
]


@parameterized_class(test_params)
class Test_warpDeriv_undeformed(unittest.TestCase):
    """AD-vs-primal-FD gate for the warping derivatives at an undeformed mesh.

    Unlike ``test_USMesh.py`` this has no reference database: the primal is
    the reference, and the test is self-verifying.
    """

    N_PROCS = 1

    def setUp(self):
        try:
            from idwarp import libidwarp  # noqa: F401
        except ImportError as err:
            raise unittest.SkipTest("Skipping because you do not have real idwarp compiled") from err

        # --- Same defaults as test_USMesh.py, so this runs the production code path ---
        self.defOpts = {
            "gridFile": None,
            "fileType": "CGNS",
            "aExp": 3.0,
            "bExp": 5.0,
            "LdefFact": 1.0,
            "alpha": 0.25,
            "errTol": 0.0005,
            "evalMode": "fast",
            "symmTol": 1e-6,
            "useRotations": True,
            "bucketSize": 8,
        }

    def test_omesh(self):
        meshOptions = copy.deepcopy(self.defOpts)
        meshOptions.update({"gridFile": os.path.join(baseDir, "../input_files/o_mesh.cgns")})
        check_warpDeriv_at_undeformed_mesh(self, meshOptions)

    def test_comesh(self):
        meshOptions = copy.deepcopy(self.defOpts)
        meshOptions.update({"gridFile": os.path.join(baseDir, "../input_files/co_mesh.cgns")})
        check_warpDeriv_at_undeformed_mesh(self, meshOptions)
