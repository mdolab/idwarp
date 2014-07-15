# this is a simple script to test the pyWarpUstruct module
import sys
import numpy
from mpi4py import MPI

# Create the global communicator
gcomm = MPI.COMM_WORLD



from pywarpustruct import USMesh

meshOptions = {}

casename = 'mdo_tutorial_node_bcs.cgns'
mesh1 = USMesh(casename, comm=gcomm, meshOptions=meshOptions)

casename = 'mdo_tutorial_face_bcs.cgns'
mesh2 = USMesh(casename, comm=gcomm, meshOptions=meshOptions)
mesh2.updateGridMetrics()
