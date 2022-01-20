# =============================================================================
# Standard Python modules
# =============================================================================
import os

# =============================================================================
# External Python modules
# =============================================================================
import numpy as np
import unittest
from mpi4py import MPI
from parameterized import parameterized_class

# =============================================================================
# Extension modules
# =============================================================================
from idwarp import MultiUSMesh, MultiUSMesh_C
from baseclasses import BaseRegTest

baseDir = os.path.dirname(os.path.abspath(__file__))  # Path to current directory


def eval_warp(handler, test_name, gridFile, optionsDict, isComplex, N_PROCS):
    # Create warping object
    if isComplex:
        # Check if the complex verision of the code has been built
        try:
            from idwarp import libidwarp_cs  # noqa: F401

            h = 1e-200
            mesh = MultiUSMesh_C(gridFile, optionsDict)
        except ImportError:
            raise unittest.SkipTest("Skipping because you do not have complex idwarp compiled")
    else:
        try:
            from idwarp import libidwarp  # noqa: F401

            mesh = MultiUSMesh(gridFile, optionsDict)
        except ImportError:
            raise unittest.SkipTest("Skipping because you do not have real idwarp compiled")

    # Get mesh data from JSON file to avoid needing ADflow for the test
    jsonFile = os.path.join(baseDir, "../input_files/onera_m6.json")
    meshData = BaseRegTest.readRefJSON(jsonFile)

    # Set the mesh surface from data
    mesh.setExternalMeshIndices(meshData["meshInd"])
    mesh.setSurfaceDefinition(meshData["pts"], meshData["conn"], meshData["faceSizes"], meshData["cgnsBlockIDs"])

    # Get the sum of the initial volume coordinates
    volNodesList = mesh.getWarpGrid()[0]
    sums = sumVolCoords(volNodesList)

    # Test the initial sums
    handler.root_add_val(f"{test_name} - Sums of Initial Volume Coords:", sums, tol=1e-14)

    # Get initial surface mesh
    coords0 = mesh.getSurfaceCoordinates()

    # Displace the surface coordinates with a shearing sweep
    new_coords = coords0.copy()
    for i in range(len(coords0)):
        span = coords0[i, 1]
        new_coords[i, 0] += 0.25 * span
    mesh.setSurfaceCoordinates(new_coords)

    # Warp the volume mesh
    mesh.warpMesh()

    # Get the sum of the warped volume coordinates
    volNodesList = mesh.getWarpGrid()[0]
    sums = sumVolCoords(volNodesList)

    # Test the warped sums
    handler.root_add_val(f"{test_name} - Sums of Warped Volume Coords:", sums, tol=1e-14)

    # Create a volume seed vector to test the mesh warping derivatives
    # We divide by the number of procs because we are not partitioning the volume mesh like ADflow would
    dXv = np.linspace(0, 1.0, len(meshData["meshInd"])) / N_PROCS

    if not isComplex:
        # Compute the warping derivatives (surface seeds)
        dXs = mesh.warpDeriv(dXv)

        # Sum the derivatives
        val = MPI.COMM_WORLD.reduce(np.sum(dXs.flatten()), op=MPI.SUM)

    else:
        # Add a complex perturbation on all surface nodes simultaneously
        for i in range(len(coords0)):
            new_coords[i, :] += 1j * h

        # Set the complex surface coordinates
        mesh.setSurfaceCoordinates(new_coords)

        # Warp the volume mesh
        mesh.warpMesh()

        # Sum the derivatives
        solverGrid = mesh.getSolverGrid()
        deriv = np.imag(solverGrid) / h
        val = MPI.COMM_WORLD.reduce(np.dot(dXv, deriv), op=MPI.SUM)

    handler.root_add_val(f"{test_name} - Sum of dXs:", val, tol=1e-14)


def sumVolCoords(volNodesList):
    """
    Helper function to take the sum of the volume coordinates for each USMesh instance.
    """

    numMeshes = len(volNodesList)
    dtype = volNodesList[0].dtype
    sums = np.zeros(numMeshes, dtype=dtype)

    for i in range(numMeshes):
        sums[i] = MPI.COMM_WORLD.reduce(np.sum(volNodesList[i]), op=MPI.SUM)

    return sums


# Set up parameterized variables
test_params = [
    {"N_PROCS": 1, "isComplex": False},
    {"N_PROCS": 2, "name": "parallel", "isComplex": False},
    {"N_PROCS": 1, "name": "complex", "isComplex": True},
    {"N_PROCS": 2, "name": "parallel_complex", "isComplex": True},
]


@parameterized_class(test_params)
class Test_MultiUSMesh(unittest.TestCase):
    """
    testflo does not run the ``train_`` functions by default.
    Use ``testflo -m "train_*" -n 1`` to regenerate the .ref files.
    """

    N_PROCS = 1

    def train_onera_m6(self):
        # Skip training complex tests, parallel tests, and the extra test that runs because of a bug in parameterized
        if not hasattr(self, "isComplex") or self.isComplex or self.N_PROCS != 1:
            raise unittest.SkipTest()
        else:
            self.test_onera_m6(train=True)

    def test_onera_m6(self, train=False):
        ref_file = os.path.join(baseDir, "ref/test_onera_m6.ref")
        with BaseRegTest(ref_file, train=train) as handler:

            # Test the ONERA M6 overset mesh
            test_name = "Test_onera_m6"
            gridFile = os.path.join(baseDir, "../input_files/onera_m6.cgns")

            meshOptions = {
                "gridFile": gridFile,
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

            optionsDict = {
                "near_tip_vol_L3": meshOptions,
                "near_wing_vol_L3": meshOptions,
            }

            # Call warping test function
            eval_warp(
                handler=handler,
                test_name=test_name,
                gridFile=gridFile,
                optionsDict=optionsDict,
                isComplex=self.isComplex,
                N_PROCS=self.N_PROCS,
            )


if __name__ == "__main__":
    unittest.main()
