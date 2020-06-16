# this is a simple script to test the IDWarp module
import sys
import numpy
from mpi4py import MPI
from idwarp import USMesh

options = {
    'gridFile':'../../input_files/mdo_tutorial_face_bcs.cgns',
    'fileType':'CGNS',
    'symmetryPlanes':[[[0,0,0], [0,0,1]]],
    'aExp': 3.0,
    'bExp': 5.0,
    'LdefFact':1.0,
    'alpha':0.25,
    'errTol':0.0005,
    'evalMode':'fast',
    'symmTol':1e-6,
    'useRotations':True,
    'bucketSize':8,
}

# Create the mesh object
mesh = USMesh(options=options, comm=MPI.COMM_WORLD)

coords0 = mesh.getSurfaceCoordinates()

new_coords = coords0.copy()
for i in range(len(coords0)):
    span = coords0[i,2]
    new_coords[i,0] += .05*span
    new_coords[i,1] += .05*span
    new_coords[i,2] += .15*span

# Reset the newly computed surface coordiantes
mesh.setSurfaceCoordinates(new_coords)

# Actually run the mesh warping
mesh.warpMesh()

# Write the new grid file.
mesh.writeGrid('warped.cgns')
