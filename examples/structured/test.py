# this is a simple script to test the IDWarp module
from mpi4py import MPI
from idwarp import USMesh

options = {
    "fileType": "CGNS",
    "gridFile": "../../input_files/o_mesh.cgns",
    "aExp": 3.0,
    "bExp": 5.0,
    "LdefFact": 100.0,
    "alpha": 0.25,
    "errTol": 0.0001,
    "evalMode": "fast",
    "symmTol": 1e-6,
    "useRotations": True,
    "bucketSize": 8,
}

# Create the mesh object
mesh = USMesh(options=options, comm=MPI.COMM_WORLD)

coords0 = mesh.getSurfaceCoordinates()

new_coords = coords0.copy()
for i in range(len(coords0)):
    new_coords[i, :] *= 1.1

# Reset the newly computed surface coordiantes
mesh.setSurfaceCoordinates(new_coords)

# Actually run the mesh warping
mesh.warpMesh()

# Write the new grid file.
mesh.writeGrid("warped.cgns")
