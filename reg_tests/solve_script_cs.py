############################################################
# DO NOT USE THIS SCRIPT AS A REFERENCE FOR HOW TO USE IDWarp
# THIS SCRIPT USES PRIVATE INTERNAL FUNCTIONALITY THAT IS
# SUBJECT TO CHANGE!!
############################################################

# ======================================================================
#         Imports
# ======================================================================
import sys
import os
import copy
import numpy

# from petsc4py import PETSc
from mdo_regression_helper import reg_write

if "complex" in sys.argv:
    from idwarp import USMesh_C as USMesh
else:
    from idwarp import USMesh

from mpi4py import MPI

# First thing we will do is define a complete set of default options
# that will be reused as we do differnt tests.  These are the default
# options as July 15, 2015.

defOpts = {
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


def printHeader(testName):
    if MPI.COMM_WORLD.rank == 0:
        print("+" + "-" * 78 + "+")
        print("| Test Name: " + "%-66s" % testName + "|")
        print("+" + "-" * 78 + "+")


h = 1e-40

# NOTE: we no longer run test1 in idwarp. test1 has been moved to pyofm and will be tested from there.
def test1():
    # Test the Ahmed body openfoam mesh
    sys.stdout.flush()
    # change directory to the correct test case
    os.chdir("../input_files/ahmedBodyMesh/")

    file_name = os.getcwd()

    meshOptions = copy.deepcopy(defOpts)

    meshOptions.update(
        {
            "gridFile": file_name,
            "fileType": "openfoam",
            "symmetryPlanes": [[[0, 0, 0], [0, 1, 0]]],
        }
    )
    # Create warping object
    mesh = USMesh(options=meshOptions, debug=True)

    # Extract Surface Coordinates
    coords0 = mesh.getSurfaceCoordinates()

    vCoords = mesh.getCommonGrid()
    val = MPI.COMM_WORLD.reduce(numpy.sum(vCoords.flatten()), op=MPI.SUM)
    if MPI.COMM_WORLD.rank == 0:
        print("Sum of vCoords Inital:")
        reg_write(val, 1e-8, 1e-8)

    new_coords = coords0.copy()
    # Do a stretch:
    for i in range(len(coords0)):
        length = coords0[i, 2]
        new_coords[i, 0] += 0.05 * length

    # Reset the newly computed surface coordiantes
    mesh.setSurfaceCoordinates(new_coords)
    mesh.warpMesh()

    vCoords = mesh.getWarpGrid()
    val = MPI.COMM_WORLD.reduce(numpy.sum(vCoords.flatten()), op=MPI.SUM)
    if MPI.COMM_WORLD.rank == 0:
        print("Sum of vCoords Warped:")
        reg_write(val, 1e-8, 1e-8)

    # Create a dXv vector to do test the mesh warping with:
    dXv_warp = numpy.linspace(0, 1.0, mesh.warp.griddata.warpmeshdof)

    if "complex" not in sys.argv:

        if MPI.COMM_WORLD.rank == 0:
            print("Computing Warp Deriv")
        mesh.warpDeriv(dXv_warp, solverVec=False)
        dXs = mesh.getdXs()

        val = MPI.COMM_WORLD.reduce(numpy.sum(dXs.flatten()), op=MPI.SUM)
        if MPI.COMM_WORLD.rank == 0:
            print("Sum of dxs:")
            reg_write(val, 1e-8, 1e-8)
    else:

        # add a complex perturbation on all surface nodes simultaneously:
        for i in range(len(coords0)):
            new_coords[i, :] += h * 1j

        # Reset the newly computed surface coordiantes
        mesh.setSurfaceCoordinates(new_coords)
        mesh.warpMesh()

        vCoords = mesh.getWarpGrid()
        deriv = numpy.imag(vCoords) / h
        deriv = numpy.dot(dXv_warp, deriv)
        val = MPI.COMM_WORLD.reduce(numpy.sum(deriv), op=MPI.SUM)

        if MPI.COMM_WORLD.rank == 0:
            print("Sum of dxs:")
            reg_write(val, 1e-8, 1e-8)

    del mesh
    # os.system('rm -fr *.cgns *.dat')
    # change back to the original directory
    os.chdir("../../reg_tests")


def test2():
    # Test the MDO tutorial co mesh
    file_name = "../input_files/co_mesh.cgns"

    meshOptions = copy.deepcopy(defOpts)

    meshOptions.update({"gridFile": file_name, "fileType": "cgns"})
    # Create warping object
    mesh = USMesh(options=meshOptions)

    # Extract Surface Coordinates
    coords0 = mesh.getSurfaceCoordinates()

    vCoords = mesh.getCommonGrid()
    val = MPI.COMM_WORLD.reduce(numpy.sum(vCoords.flatten()), op=MPI.SUM)
    if MPI.COMM_WORLD.rank == 0:
        print("Sum of vCoords Inital:")
        reg_write(val, 1e-8, 1e-8)

    new_coords = coords0.copy()
    # Do a shearing sweep deflection:
    for i in range(len(coords0)):
        span = coords0[i, 2]
        new_coords[i, 0] += 0.05 * span

    # Reset the newly computed surface coordiantes
    mesh.setSurfaceCoordinates(new_coords)
    mesh.warpMesh()

    # Get the sum of the warped coordinates
    # vCoords = mesh.getSolverGrid()
    vCoords = mesh.getWarpGrid()

    val = MPI.COMM_WORLD.reduce(numpy.sum(vCoords.flatten()), op=MPI.SUM)
    if MPI.COMM_WORLD.rank == 0:
        print("Sum of vCoords Warped:")
        reg_write(val, 1e-5, 1e-5)

    # Create a dXv vector to do test the mesh warping with:
    dXv_warp = numpy.linspace(0, 1.0, mesh.warp.griddata.warpmeshdof)

    if "complex" not in sys.argv:
        if MPI.COMM_WORLD.rank == 0:
            print("Computing Warp Deriv")
        mesh.warpDeriv(dXv_warp, solverVec=False)
        dXs = mesh.getdXs()

        val = MPI.COMM_WORLD.reduce(numpy.sum(dXs.flatten()), op=MPI.SUM)
        if MPI.COMM_WORLD.rank == 0:
            print("Sum of dxs:")
            reg_write(val, 1e-5, 1e-5)
    else:

        # add a complex perturbation on all surface nodes simultaneously:
        for i in range(len(coords0)):
            new_coords[i, :] += h * 1j

        # Reset the newly computed surface coordiantes
        mesh.setSurfaceCoordinates(new_coords)
        mesh.warpMesh()

        vCoords = mesh.getWarpGrid()
        deriv = numpy.imag(vCoords) / h
        deriv = numpy.dot(dXv_warp, deriv)
        val = MPI.COMM_WORLD.reduce(numpy.sum(deriv), op=MPI.SUM)

        if MPI.COMM_WORLD.rank == 0:
            print("Sum of dxs:")
            reg_write(val, 1e-8, 1e-8)

    del mesh
    # os.system('rm -fr *.cgns *.dat')


def test3():
    # Test the mdo tutorial o mesh
    file_name = "../input_files/o_mesh.cgns"

    meshOptions = copy.deepcopy(defOpts)

    meshOptions.update({"gridFile": file_name, "fileType": "cgns"})
    # Create warping object
    mesh = USMesh(options=meshOptions)

    # Extract Surface Coordinates
    coords0 = mesh.getSurfaceCoordinates()

    vCoords = mesh.getCommonGrid()
    val = MPI.COMM_WORLD.reduce(numpy.sum(vCoords.flatten()), op=MPI.SUM)
    if MPI.COMM_WORLD.rank == 0:
        print("Sum of vCoords Inital:")
        reg_write(val, 1e-8, 1e-8)

    new_coords = coords0.copy()
    # Do a shearing sweep deflection:
    for i in range(len(coords0)):
        span = coords0[i, 2]
        new_coords[i, 0] += 0.05 * span

    # Reset the newly computed surface coordiantes
    mesh.setSurfaceCoordinates(new_coords)
    mesh.warpMesh()

    # Get the sum of the warped coordinates
    # vCoords = mesh.getSolverGrid()
    vCoords = mesh.getWarpGrid()

    val = MPI.COMM_WORLD.reduce(numpy.sum(vCoords.flatten()), op=MPI.SUM)
    if MPI.COMM_WORLD.rank == 0:
        print("Sum of vCoords Warped:")
        reg_write(val, 1e-8, 1e-8)

    # Create a dXv vector to do test the mesh warping with:
    dXv_warp = numpy.linspace(0, 1.0, mesh.warp.griddata.warpmeshdof)

    if "complex" not in sys.argv:
        if MPI.COMM_WORLD.rank == 0:
            print("Computing Warp Deriv")
        mesh.warpDeriv(dXv_warp, solverVec=False)
        dXs = mesh.getdXs()

        val = MPI.COMM_WORLD.reduce(numpy.sum(dXs.flatten()), op=MPI.SUM)
        if MPI.COMM_WORLD.rank == 0:
            print("Sum of dxs:")
            reg_write(val, 1e-8, 1e-8)
    else:

        # add a complex perturbation on all surface nodes simultaneously:
        for i in range(len(coords0)):
            new_coords[i, :] += h * 1j

        # Reset the newly computed surface coordiantes
        mesh.setSurfaceCoordinates(new_coords)
        mesh.warpMesh()

        vCoords = mesh.getWarpGrid()
        deriv = numpy.imag(vCoords) / h
        deriv = numpy.dot(dXv_warp, deriv)
        val = MPI.COMM_WORLD.reduce(numpy.sum(deriv), op=MPI.SUM)

        if MPI.COMM_WORLD.rank == 0:
            print("Sum of dxs:")
            reg_write(val, 1e-8, 1e-8)


if __name__ == "__main__":
    if len(sys.argv) == 1 or (len(sys.argv) == 2 and "complex" in sys.argv):
        # NOTE: commented out test1. test1 has been moved to pyofom and will be tested from there.
        # test1()
        test2()
        test3()

    else:
        # Run individual ones
        # NOTE: commented out test1. test1 has been moved to pyofom and will be tested from there.
        # if 'test1' in sys.argv:
        #    test1()
        if "test2" in sys.argv:
            test2()
        if "test3" in sys.argv:
            test3()
