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

baseDir = os.path.dirname(os.path.abspath(__file__))


def eval_warp(handler, test_name, meshOptions):
    # Create warping object
    mesh = USMesh(options=meshOptions)

    # Extract Surface Coordinates
    coords0 = mesh.getSurfaceCoordinates()

    vCoords = mesh.getCommonGrid()
    val = MPI.COMM_WORLD.reduce(numpy.sum(vCoords.flatten()), op=MPI.SUM)
    handler.root_add_val(f"{test_name} - Sum of vCoords Inital:", val, tol=1e-8)

    new_coords = coords0.copy()
    # Do a shearing sweep deflection:
    if test_name == "Test_inflate_cube":
        for i in range(len(coords0)):
            new_coords[i, 0] *= 1.1
            new_coords[i, 1] *= 1.2
            new_coords[i, 1] *= 1.3
    else:
        for i in range(len(coords0)):
            span = coords0[i, 2]
            new_coords[i, 0] += 0.05 * span

    # Reset the newly computed surface coordiantes
    mesh.setSurfaceCoordinates(new_coords)
    mesh.warpMesh()

    # Get the sum of the warped coordinates
    vCoords = mesh.getWarpGrid()

    val = MPI.COMM_WORLD.reduce(numpy.sum(vCoords.flatten()), op=MPI.SUM)
    handler.root_add_val(f"{test_name} - Sum of vCoords Warped:", val, tol=1e-8)

    # Create a dXv vector to do test the mesh warping with:
    dXv_warp = numpy.linspace(0, 1.0, mesh.warp.griddata.warpmeshdof)

    # Computing Warp Derivatives
    mesh.warpDeriv(dXv_warp, solverVec=False)
    dXs = mesh.getdXs()

    val = MPI.COMM_WORLD.reduce(numpy.sum(dXs.flatten()), op=MPI.SUM)
    handler.root_add_val(f"{test_name} - Sum of dxs:", val, tol=1e-8)

    # Verifying Warp Deriv
    # TODO: Do we want to keep / exploit this automatic deriv check?
    mesh.verifyWarpDeriv(dXv_warp, solverVec=False, dofStart=0, dofEnd=5)


test_params = [{"N_PROCS": 1}, {"N_PROCS": 2, "name": "parallel"}]


@parameterized_class(test_params)
class Test_USmesh(unittest.TestCase):

    N_PROCS = 1

    def setUp(self):
        # TODO keep this explicit set of options? Or use the default ones?
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
        self.ref_app = ""

        if self.N_PROCS > 1:
            self.ref_app = "_par"
            print(self.ref_app)

    def train_comesh(self, train=True):
        self.test_comesh(train=train)

    def test_comesh(self, train=False):
        print(self.ref_app)
        ref_file = os.path.join(baseDir, f"ref/test_comesh{self.ref_app}.ref")
        with BaseRegTest(ref_file, train=train) as handler:
            # Test the mdo tutorial co mesh
            test_name = "Test_co_mesh"
            file_name = os.path.join(baseDir, "../input_files/co_mesh.cgns")

            meshOptions = copy.deepcopy(self.defOpts)

            meshOptions.update({"gridFile": file_name, "fileType": "cgns"})

            eval_warp(handler=handler, test_name=test_name, meshOptions=meshOptions)

    def train_omesh(self, train=True):
        self.test_omesh(train=train)

    def test_omesh(self, train=False):
        ref_file = os.path.join(baseDir, f"ref/test_omesh{self.ref_app}.ref")
        with BaseRegTest(ref_file, train=train) as handler:
            # Test the mdo tutorial o mesh
            test_name = "Test_o_mesh"
            file_name = os.path.join(baseDir, "../input_files/o_mesh.cgns")

            meshOptions = copy.deepcopy(self.defOpts)

            meshOptions.update({"gridFile": file_name, "fileType": "cgns"})

            eval_warp(handler=handler, test_name=test_name, meshOptions=meshOptions)

    def train_sym_mesh(self, train=True):
        self.test_sym_mesh(train=train)

    def test_sym_mesh(self, train=False):
        ref_file = os.path.join(baseDir, f"ref/test_sym_mesh{self.ref_app}.ref")
        with BaseRegTest(ref_file, train=train) as handler:
            # Test the MDO tutorial h mesh
            test_name = "Test_h_mesh"
            file_name = os.path.join(baseDir, "../input_files/mdo_tutorial_face_bcs.cgns")

            meshOptions = copy.deepcopy(self.defOpts)

            meshOptions.update(
                {
                    "gridFile": file_name,
                    "fileType": "cgns",
                    "symmetryPlanes": [[[0, 0, 0], [0, 0, 1]]],
                }
            )

            eval_warp(handler=handler, test_name=test_name, meshOptions=meshOptions)

    def train_inflate_cube(self, train=True):
        self.test_inflate_cube(train=train)

    def test_inflate_cube(self, train=False):
        ref_file = os.path.join(baseDir, f"ref/test_inflate_cube{self.ref_app}.ref")
        with BaseRegTest(ref_file, train=train) as handler:
            # Test the MDO tutorial h mesh
            test_name = "Test_inflate_cube"
            file_name = os.path.join(baseDir, "../input_files/symm_block.cgns")

            meshOptions = copy.deepcopy(self.defOpts)

            meshOptions.update(
                {
                    "gridFile": file_name,
                    "fileType": "cgns",
                }
            )

            eval_warp(handler=handler, test_name=test_name, meshOptions=meshOptions)
