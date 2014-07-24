# this is a simple script to test the pyWarpUstruct module
import sys
import numpy
from mpi4py import MPI

# Create the global communicator
gcomm = MPI.COMM_WORLD

from pywarpustruct import USMesh

meshOptions = {}

# casename = 'mdo_tutorial_node_bcs.cgns'
# mesh1 = USMesh(casename, comm=gcomm, meshOptions=meshOptions)

casename = 'mdo_tutorial_face_bcs.cgns'
mesh = USMesh(casename, comm=gcomm, meshOptions=meshOptions)
print 'mesh initialized',mesh.familyGroup
# First see what families are in the file:
mesh.printFamilyList()

# Add two family groups: 
mesh.addFamilyGroup('upper_surface',['wing_up']) 
mesh.addFamilyGroup('lower_surface',['wing_low'])
mesh.addFamilyGroup('full_surface',['wing_low','wing_up','tip_up','tip_low'])
print 'families added',mesh.familyGroup

coords0 = mesh.getSurfaceCoordinates('full_surface')
print 'updating coordinates'
new_coords = coords0.copy()
for i in xrange(len(coords0)):
    span = coords0[i,2]
    new_coords[i,0] += .05*span
# end for
print 'setting new coords'
# Reset the newly computed surface coordiantes
mesh.setSurfaceCoordinates('full_surface',new_coords)

mesh.warpMesh()
sys.exit(0)

# mesh2.preprocessSurfaces()
# mesh2.updateGridMetrics()
