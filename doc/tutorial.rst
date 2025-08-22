.. _tutorial:

Tutorial
========

This is a short tutorial explaining how to get started with
IDWarp. We will start by presenting a simple stand-alone
warping script and explaining the various parts. ::

  from mpi4py import MPI
  from idwarp import USMesh

  options = {
    'gridFile':'../../input_files/o_mesh.cgns',
    'fileType':'cgns',
    'specifiedSurfaces':None,
    'symmetrySurfaces':None,
    'symmetryPlanes':[],
    'aExp': 3.0,
    'bExp': 5.0,
    'LdefFact':1.0,
    'alpha':0.25,
    'errTol':0.0001,
    'evalMode':'fast',
    'useRotations':True,
    'zeroCornerRotations':True,
    'cornerAngle':30.0,
    'bucketSize':8,
  }

  # Create the mesh object
  mesh = USMesh(options=options, comm=MPI.COMM_WORLD)

  # Extract all coordinates
  coords0 = mesh.getSurfaceCoordinates()

  # Modify the coordinates as required
  newCoords = coords0.copy()
  for i in range(len(coords0)):
      newCoords[i, :] *= 1.1

  # Reset the newly computed surface coordinates
  mesh.setSurfaceCoordinates(newCoords)

  # Actually run the mesh warping
  mesh.warpMesh()

  # Write the new grid file.
  mesh.writeGrid('warped.cgns')

The first two lines of code imports mpi4py and the mesh warping module ::

  from mpi4py import MPI
  from idwarp import USMesh

The next chunk lists the options for warping. The options are
explained in :ref:`idwarp_options`. ::

  options = {
    'gridFile':'../../input_files/o_mesh.cgns',
    'fileType':'cgns',
    'specifiedSurfaces':None,
    'symmetrySurfaces':None,
    'symmetryPlanes':[],
    'aExp': 3.0,
    'bExp': 5.0,
    'LdefFact':1.0,
    'alpha':0.25,
    'errTol':0.0001,
    'evalMode':'fast',
    'useRotations':True,
    'zeroCornerRotations':True,
    'cornerAngle':30.0,
    'bucketSize':8,
  }

Next we create the actual mesh object itself ::

 # Create the mesh object
 mesh = USMesh(options=options, comm=MPI.COMM_WORLD)

Note that we have explicitly passed in the MPI intracommunicator on
which we want to create the object. If the 'comm' keyword argument is
not given, it will default to MPI.COMM_WORLD. Therefor this example,
mpi4py is not strictly required to be imported in the run script.

Next we request the surface coordinates from the mesh object. These
will correspond to coordinates in the 'specifiedSurfaces' option. ::

  # Extract all coordinates
  coords0 = mesh.getSurfaceCoordinates()

coords0 is a numpy array of size (N,3). It is now up to the user to
manipulate these coordinates however they wish for this example we
simply loop over all coordinates and uniformly scale by a factor of 1.1::

  new_coords = coords0.copy()
  for i in range(len(coords0)):
      new_coords[i, :] *= 1.1

Once the new set of coordinates have been determined, return them to
the mesh warping object with the following command. ::

  # Reset the newly computed surface coordiantes
  mesh.setSurfaceCoordinates(new_coords)

Note that the shape of 'new_coords' must be identical to the coords0
array that was originally provided by the warping. Next we run the
actual mesh warp using ::

  # Actually run the mesh warping
  mesh.warpMesh()

And finally to produce an updated grid file we can write the grid::

  # Write the new grid file.
  mesh.writeGrid('warped.cgns')

The warped grid file 'warped.cgns' will contain all the boundary
condition/connectivity/auxiliary information as the original cgns
file. Only the coordinates are updated to their new positions.
