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


class Test(unittest.TestCase):

    N_PROCS = 1

    def setUp(self):
        self.ref_file = os.path.join(baseDir, "ref/test_.ref")

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

    def train_comesh(self, train=True):
        self.test_comesh(train=train)
        # handler.writeRef()

    def test_comesh(self, train=False):
        with BaseRegTest(self.ref_file, train=train) as handler:
            test_name = "Test_co_mesh"
            file_name = os.path.join(baseDir, "../input_files/co_mesh.cgns")

            meshOptions = copy.deepcopy(self.defOpts)

            meshOptions.update({"gridFile": file_name, "fileType": "cgns"})

            eval_warp(handler=handler, test_name=test_name, meshOptions=meshOptions)

            if train is True:
                handler.writeRef()
