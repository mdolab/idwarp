# this is a simple script to test the pyWarpUstruct module
import sys
import numpy
from mpi4py import MPI
from pywarpustruct import USMesh

options = {
    'gridFile':'../input_files/co_mesh.cgns',
    'fileType':'CGNS',
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

# First see what families are in the file:
mesh.printFamilyList()

# Add two family groups: 
mesh.addFamilyGroup('upper_surface',['wing_up']) 
mesh.addFamilyGroup('lower_surface',['wing_low'])
mesh.addFamilyGroup('full_surface',['wing_low','wing_up'])

coords0 = mesh.getSurfaceCoordinates('full_surface')

new_coords = coords0.copy()
for i in xrange(len(coords0)):
    span = coords0[i,2]
    new_coords[i,0] += .05*span
    new_coords[i,1] += .05*span
    new_coords[i,2] += .15*span

# Reset the newly computed surface coordiantes
mesh.setSurfaceCoordinates(new_coords, 'full_surface')

# Actually run the mesh warping
mesh.warpMesh()

# Write the new grid file.
mesh.writeGrid('warped.cgns')
