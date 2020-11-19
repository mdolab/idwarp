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
from baseclasses import BaseRegTest

baseDir = os.path.dirname(os.path.abspath(__file__))  # Path to current folder


def eval_warp(handler, test_name, meshOptions):
    # --- Create warping object ---
    mesh = USMesh(options=meshOptions)

    # --- Extract Surface Coordinates ---
    coords0 = mesh.getSurfaceCoordinates()

    vCoords = mesh.getCommonGrid()
    val = MPI.COMM_WORLD.reduce(numpy.sum(vCoords.flatten()), op=MPI.SUM)
    handler.root_add_val(f"{test_name} - Sum of vCoords Inital:", val, tol=1e-8)

    new_coords = coords0.copy()

    # --- Setting the coordinates displacement for surface mesh ---
    if test_name == "Test_inflate_cube":
        # This specific test displaces the surface in every direction
        for i in range(len(coords0)):
            new_coords[i, 0] *= 1.1
            new_coords[i, 1] *= 1.2
            new_coords[i, 1] *= 1.3  # TODO this should be [i,2], right?
    else:
        # Do a shearing sweep deflection:
        for i in range(len(coords0)):
            span = coords0[i, 2]
            new_coords[i, 0] += 0.05 * span

    # --- Reset the newly computed surface coordiantes
    mesh.setSurfaceCoordinates(new_coords)
    mesh.warpMesh()

    # --- Get the sum of the warped coordinates ---
    vCoords = mesh.getWarpGrid()

    val = MPI.COMM_WORLD.reduce(numpy.sum(vCoords.flatten()), op=MPI.SUM)
    handler.root_add_val(f"{test_name} - Sum of vCoords Warped:", val, tol=1e-8)

    # --- Create a dXv vector to do test the mesh warping with ---
    dXv_warp = numpy.linspace(0, 1.0, mesh.warp.griddata.warpmeshdof)

    # --- Computing Warp Derivatives ---
    mesh.warpDeriv(dXv_warp, solverVec=False)
    dXs = mesh.getdXs()

    val = MPI.COMM_WORLD.reduce(numpy.sum(dXs.flatten()), op=MPI.SUM)
    handler.root_add_val(f"{test_name} - Sum of dxs:", val, tol=1e-8)

    # --- Verifying Warp Deriv ---
    # TODO: Do we want to keep / exploit this automatic deriv check?
    mesh.verifyWarpDeriv(dXv_warp, solverVec=False, dofStart=0, dofEnd=5)


# --- Set up variables to be overridden at every parametrized_class loop ---
test_params = [{"N_PROCS": 1}, {"N_PROCS": 2, "name": "parallel"}]


@parameterized_class(test_params)
class Test_USmesh(unittest.TestCase):
    """
    Defines the tests for IDwarp.
    'train_' methods are not run by default by testflo.
    Use `$ testflo ./ -m "train_*"` to re-generate the .ref files used as benchmark.
    Note that, in the current state, you will get a few FAIL messages when training - that's due to an issue with parameterize, the .ref files will be generated correctly.

    Also note that the `handler` BaseRegTest class will automatically take care to benchmark or train the tests (including file IO)
    """

    N_PROCS = 1

    def setUp(self):
        # TODO keep this explicit set of options? Or use the default ones?
        # --- Reference options to be used by every test ---
        self.defOpts = {
            "gridFile": None,
            "aExp": 3.0,
            "bExp": 5.0,
            "LdefFact": 1.0,
            "alpha": 0.25,
            "errTol": 0.0005,
            "evalMode": "fast",
            "symmTol": 1e-6,
            "useRotations": True,
            "bucketSize": 8,
            "fileType": None,
        }
        # --- Setting the ref file for parallel tests ---
        self.ref_app = ""

        if self.N_PROCS > 1:
            # Appendix for the parallel reference files - appended to every ref_file instance
            self.ref_app = "_par"

    def train_comesh(self, train=True):
        try:
            self.test_comesh(train=train)
        except AttributeError:
            raise unittest.SkipTest()

    def test_comesh(self, train=False):
        ref_file = os.path.join(baseDir, f"ref/test_comesh{self.ref_app}.ref")
        with BaseRegTest(ref_file, train=train) as handler:
            # Test the C-mesh in the input files
            test_name = "Test_co_mesh"
            file_name = os.path.join(baseDir, "../input_files/co_mesh.cgns")

            meshOptions = copy.deepcopy(self.defOpts)

            meshOptions.update({"gridFile": file_name, "fileType": "cgns"})

            # --- Call shared computation --
            eval_warp(handler=handler, test_name=test_name, meshOptions=meshOptions)

    def train_omesh(self, train=True):
        try:
            self.test_omesh(train=train)
        except AttributeError:
            raise unittest.SkipTest()

    def test_omesh(self, train=False):
        ref_file = os.path.join(baseDir, f"ref/test_omesh{self.ref_app}.ref")
        with BaseRegTest(ref_file, train=train) as handler:
            # --- Test the O-mesh ---
            test_name = "Test_o_mesh"
            file_name = os.path.join(baseDir, "../input_files/o_mesh.cgns")

            meshOptions = copy.deepcopy(self.defOpts)

            meshOptions.update({"gridFile": file_name, "fileType": "cgns"})

            # --- Call shared computation ---
            eval_warp(handler=handler, test_name=test_name, meshOptions=meshOptions)

    def train_sym_mesh(self, train=True):
        try:
            self.test_sym_mesh(train=train)
        except AttributeError:
            raise unittest.SkipTest()

    def test_sym_mesh(self, train=False):
        ref_file = os.path.join(baseDir, f"ref/test_sym_mesh{self.ref_app}.ref")
        with BaseRegTest(ref_file, train=train) as handler:
            # Test the symmetric mesh
            test_name = "Test_sym_mesh"
            file_name = os.path.join(baseDir, "../input_files/mdo_tutorial_face_bcs.cgns")

            meshOptions = copy.deepcopy(self.defOpts)

            meshOptions.update(
                {
                    "gridFile": file_name,
                    "fileType": "cgns",
                    "symmetryPlanes": [[[0, 0, 0], [0, 0, 1]]],
                }
            )

            # --- Call shared computation --
            eval_warp(handler=handler, test_name=test_name, meshOptions=meshOptions)

    def train_inflate_cube(self, train=True):
        try:
            self.test_inflate_cube(train=train)
        except AttributeError:
            raise unittest.SkipTest()

    def test_inflate_cube(self, train=False):
        ref_file = os.path.join(baseDir, f"ref/test_inflate_cube{self.ref_app}.ref")
        with BaseRegTest(ref_file, train=train) as handler:
            # Test the "cube" mesh
            test_name = "Test_inflate_cube"
            file_name = os.path.join(baseDir, "../input_files/symm_block.cgns")

            meshOptions = copy.deepcopy(self.defOpts)

            meshOptions.update({"gridFile": file_name, "fileType": "cgns"})

            # --- Call shared computation --
            eval_warp(handler=handler, test_name=test_name, meshOptions=meshOptions)
