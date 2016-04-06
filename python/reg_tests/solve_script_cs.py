
############################################################
# DO NOT USE THIS SCRIPT AS A REFERENCE FOR HOW TO USE pywarpustruct
# THIS SCRIPT USES PRIVATE INTERNAL FUNCTIONALITY THAT IS
# SUBJECT TO CHANGE!!
############################################################

# ======================================================================
#         Imports
# ======================================================================
from __future__ import print_function
import sys,os,copy
import numpy
from mdo_regression_helper import *
if 'complex' in sys.argv:
    from pywarpustruct import USMesh_C as USMesh
else:
    from pywarpustruct import USMesh

from mpi4py import MPI

# First thing we will do is define a complete set of default options
# that will be reused as we do differnt tests.  These are the default
# options as July 15, 2015.

defOpts = {
    'gridFile':None,
    'aExp': 3.0,
    'bExp': 5.0,
    'LdefFact':1.0,
    'alpha':0.25,
    'errTol':0.0005,
    'evalMode':'fast',
    'symmTol':1e-6,
    'useRotations':True,
    'bucketSize':8,
    'fileType':None
}


def printHeader(testName):
    if MPI.COMM_WORLD.rank == 0:
        print('+' + '-'*78 + '+')
        print('| Test Name: ' + '%-66s'%testName + '|')
        print('+' + '-'*78 + '+')

h = 1e-40

def test1():
    # Test the Ahmed body openfoam mesh
    sys.stdout.flush()
    #change directory to the correct test case
    os.chdir('../input_files/ahmedBodyMesh/')

    file_name = os.getcwd()

    meshOptions = copy.deepcopy(defOpts)

    meshOptions.update(
        {'gridFile':file_name,
         'fileType':'openfoam'
     }
    )
    # Create warping object
    mesh = USMesh(options=meshOptions)

    # Add two family groups: 
    mesh.addFamilyGroup('body', ['wall']) 
    mesh.addFamilyGroup('ground', ['lowerwall'])

    # Extract Surface Coordinates
    coords0 = mesh.getSurfaceCoordinates('body')

    vCoords = mesh.getCommonGrid()
    val = MPI.COMM_WORLD.reduce(numpy.sum(vCoords.flatten()),op=MPI.SUM)
    if MPI.COMM_WORLD.rank == 0:
        print('Sum of vCoords Inital:')
        reg_write(val,1e-8,1e-8)

    new_coords = coords0.copy()
    # Do a stretch:
    for i in xrange(len(coords0)):
        length = coords0[i,2]
        new_coords[i,0] += .05*length

    # Reset the newly computed surface coordiantes
    mesh.setSurfaceCoordinates(new_coords, 'body')
    mesh.warpMesh()

    vCoords = mesh.getWarpGrid()
    val = MPI.COMM_WORLD.reduce(numpy.sum(vCoords.flatten()),op=MPI.SUM)
    if MPI.COMM_WORLD.rank == 0:
        print('Sum of vCoords Warped:')
        reg_write(val,1e-8,1e-8)

    # Create a dXv vector to do test the mesh warping with:
    dXv_warp = numpy.linspace(0,1.0, mesh.warp.griddata.warpmeshdof)

    if not 'complex' in sys.argv:
       
        if MPI.COMM_WORLD.rank == 0:
            print('Computing Warp Deriv')
        mesh.warpDeriv(dXv_warp,solverVec=False)
        dXs = mesh.getdXs('body')
            
        val = MPI.COMM_WORLD.reduce(numpy.sum(dXs.flatten()),op=MPI.SUM)
        if MPI.COMM_WORLD.rank == 0:
            print('Sum of dxs:')
            reg_write(val,1e-8,1e-8)
    else:
               
        # add a complex perturbation on all surface nodes simultaneously:
        for i in xrange(len(coords0)):
            new_coords[i,:] += h*1j

        # Reset the newly computed surface coordiantes
        mesh.setSurfaceCoordinates(new_coords, 'body')
        mesh.warpMesh()  

        vCoords = mesh.getWarpGrid()
        deriv = numpy.imag(vCoords)/h
        deriv = numpy.dot(dXv_warp,deriv)
        val = MPI.COMM_WORLD.reduce(numpy.sum(deriv),op=MPI.SUM)
        val/=len(coords0)*3
        if MPI.COMM_WORLD.rank == 0:
            print('Sum of dxs:')
            reg_write(val,1e-8,1e-8)
        

    del mesh
    #os.system('rm -fr *.cgns *.dat')
    #change back to the original directory
    os.chdir('../../reg_tests')


def test2():
    # Test the MDO tutorial co mesh
    file_name = '../input_files/co_mesh.cgns'

    meshOptions = copy.deepcopy(defOpts)

    meshOptions.update(
        {'gridFile':file_name,
         'fileType':'cgns'
     }
    )
    # Create warping object
    mesh = USMesh(options=meshOptions)

    # Add family group
    mesh.addFamilyGroup('full_surface')
    
    # Extract Surface Coordinates
    coords0 = mesh.getSurfaceCoordinates('full_surface')

    vCoords = mesh.getCommonGrid()
    val = MPI.COMM_WORLD.reduce(numpy.sum(vCoords.flatten()),op=MPI.SUM)
    if MPI.COMM_WORLD.rank == 0:
        print('Sum of vCoords Inital:')
        reg_write(val,1e-8,1e-8)

    new_coords = coords0.copy()
    # Do a shearing sweep deflection:
    for i in xrange(len(coords0)):
        span = coords0[i,2]
        new_coords[i,0] += .05*span

    # Reset the newly computed surface coordiantes
    print('setting surface')
    mesh.setSurfaceCoordinates(new_coords, 'full_surface')
    print('warping mesh')
    mesh.warpMesh()
    
    # Get the sum of the warped coordinates
    #vCoords = mesh.getSolverGrid()
    vCoords = mesh.getWarpGrid()

    val = MPI.COMM_WORLD.reduce(numpy.sum(vCoords.flatten()),op=MPI.SUM)
    if MPI.COMM_WORLD.rank == 0:
        print('Sum of vCoords Warped:')
        reg_write(val,1e-8,1e-8)

    # Create a dXv vector to do test the mesh warping with:
    dXv_warp = numpy.linspace(0,1.0, mesh.warp.griddata.warpmeshdof)

    if not 'complex' in sys.argv:
        if MPI.COMM_WORLD.rank == 0:
            print('Computing Warp Deriv')
        mesh.warpDeriv(dXv_warp,solverVec=False)
        dXs = mesh.getdXs('full_surface')
            
        val = MPI.COMM_WORLD.reduce(numpy.sum(dXs.flatten()),op=MPI.SUM)
        if MPI.COMM_WORLD.rank == 0:
            print('Sum of dxs:')
            reg_write(val,1e-8,1e-8)
    else:
        
        # add a complex perturbation on all surface nodes simultaneously:
        for i in xrange(len(coords0)):
            new_coords[i,:] += h*1j

        # Reset the newly computed surface coordiantes
        mesh.setSurfaceCoordinates(new_coords, 'full_surface')
        mesh.warpMesh()  

        vCoords = mesh.getWarpGrid()
        deriv = numpy.imag(vCoords)/h
        deriv = numpy.dot(dXv_warp,deriv)
        val = MPI.COMM_WORLD.reduce(numpy.sum(deriv),op=MPI.SUM)
        val/=len(coords0)*3
        if MPI.COMM_WORLD.rank == 0:
            print('Sum of dxs:')
            reg_write(val,1e-8,1e-8)

    del mesh
    #os.system('rm -fr *.cgns *.dat')

def test3():
    # Test the mdo tutorial o mesh
    file_name = '../input_files/o_mesh.cgns'

    meshOptions = copy.deepcopy(defOpts)

    meshOptions.update(
        {'gridFile':file_name,
         'fileType':'cgns'
     }
    )
    # Create warping object
    mesh = USMesh(options=meshOptions)

    # Add family group
    mesh.addFamilyGroup('full_surface')
    
    # Extract Surface Coordinates
    coords0 = mesh.getSurfaceCoordinates('full_surface')

    vCoords = mesh.getCommonGrid()
    val = MPI.COMM_WORLD.reduce(numpy.sum(vCoords.flatten()),op=MPI.SUM)
    if MPI.COMM_WORLD.rank == 0:
        print('Sum of vCoords Inital:')
        reg_write(val,1e-8,1e-8)

    new_coords = coords0.copy()
    # Do a shearing sweep deflection:
    for i in xrange(len(coords0)):
        span = coords0[i,2]
        new_coords[i,0] += .05*span

    # Reset the newly computed surface coordiantes
    print('setting surface')
    mesh.setSurfaceCoordinates(new_coords, 'full_surface')
    print('warping mesh')
    mesh.warpMesh()
    
    # Get the sum of the warped coordinates
    #vCoords = mesh.getSolverGrid()
    vCoords = mesh.getWarpGrid()

    val = MPI.COMM_WORLD.reduce(numpy.sum(vCoords.flatten()),op=MPI.SUM)
    if MPI.COMM_WORLD.rank == 0:
        print('Sum of vCoords Warped:')
        reg_write(val,1e-8,1e-8)

    # Create a dXv vector to do test the mesh warping with:
    dXv_warp = numpy.linspace(0,1.0, mesh.warp.griddata.warpmeshdof)

    if not 'complex' in sys.argv:
        if MPI.COMM_WORLD.rank == 0:
            print('Computing Warp Deriv')
        mesh.warpDeriv(dXv_warp,solverVec=False)
        dXs = mesh.getdXs('full_surface')
            
        val = MPI.COMM_WORLD.reduce(numpy.sum(dXs.flatten()),op=MPI.SUM)
        if MPI.COMM_WORLD.rank == 0:
            print('Sum of dxs:')
            reg_write(val,1e-8,1e-8)
    else:
        
        # add a complex perturbation on all surface nodes simultaneously:
        for i in xrange(len(coords0)):
            new_coords[i,:] += h*1j

        # Reset the newly computed surface coordiantes
        mesh.setSurfaceCoordinates(new_coords, 'full_surface')
        mesh.warpMesh()  

        vCoords = mesh.getWarpGrid()
        deriv = numpy.imag(vCoords)/h
        deriv = numpy.dot(dXv_warp,deriv)
        val = MPI.COMM_WORLD.reduce(numpy.sum(deriv),op=MPI.SUM)
        val/=len(coords0)*3
        if MPI.COMM_WORLD.rank == 0:
            print('Sum of dxs:')
            reg_write(val,1e-8,1e-8)

    del mesh
    #os.system('rm -fr *.cgns *.dat')

# def test4():
#     # Test the MDO tutorial h mesh
#     file_name = '../examples/input/h_mesh.cgns'

#     meshOptions = copy.deepcopy(defOpts)

#     meshOptions.update(
#         {'gridFile':file_name,
#          'fileType':'cgns'
#      }
#     )
#     # Create warping object
#     mesh = USMesh(options=meshOptions)

#     # Add family group
#     mesh.addFamilyGroup('full_surface')
    
#     # Extract Surface Coordinates
#     coords0 = mesh.getSurfaceCoordinates('full_surface')

#     vCoords = mesh.getCommonGrid()
#     val = MPI.COMM_WORLD.reduce(numpy.sum(vCoords.flatten()),op=MPI.SUM)
#     if MPI.COMM_WORLD.rank == 0:
#         print('Sum of vCoords Inital:')
#         reg_write(val,1e-8,1e-8)

#     new_coords = coords0.copy()
#     # Do a shearing sweep deflection:
#     for i in xrange(len(coords0)):
#         span = coords0[i,2]
#         new_coords[i,0] += .05*span

#     # Reset the newly computed surface coordiantes
#     print('setting surface')
#     mesh.setSurfaceCoordinates(new_coords, 'full_surface')
#     print('warping mesh')
#     mesh.warpMesh()
    
#     # Get the sum of the warped coordinates
#     #vCoords = mesh.getSolverGrid()
#     vCoords = mesh.getWarpGrid()

#     val = MPI.COMM_WORLD.reduce(numpy.sum(vCoords.flatten()),op=MPI.SUM)
#     if MPI.COMM_WORLD.rank == 0:
#         print('Sum of vCoords Warped:')
#         reg_write(val,1e-8,1e-8)

#     # Create a dXv vector to do test the mesh warping with:
#     dXv_warp = numpy.linspace(0,1.0, mesh.warp.griddata.warpmeshdof)
#     if MPI.COMM_WORLD.rank == 0:
#         print('Computing Warp Deriv')
#     mesh.warpDeriv(dXv_warp,solverVec=False)
#     dXs = mesh.getdXs('full_surface')

#     val = MPI.COMM_WORLD.reduce(numpy.sum(dXs.flatten()),op=MPI.SUM)
#     if MPI.COMM_WORLD.rank == 0:
#         print('Sum of dxs:')
#         reg_write(val,1e-8,1e-8)
        
#     if MPI.COMM_WORLD.rank == 0:
#         print('Verifying Warp Deriv')
#     mesh.verifyWarpDeriv(dXv_warp,solverVec=False,dofStart=0,dofEnd=5)

#     del mesh
#     #os.system('rm -fr *.cgns *.dat')

def test5():
    # Test the MDO tutorial h mesh
    file_name = '../input_files/mdo_tutorial_face_bcs.cgns'

    meshOptions = copy.deepcopy(defOpts)

    meshOptions.update(
        {'gridFile':file_name,
         'fileType':'cgns'
     }
    )
    # Create warping object
    mesh = USMesh(options=meshOptions)

    # Add family group
    mesh.addFamilyGroup('full_surface')
    
    # Extract Surface Coordinates
    coords0 = mesh.getSurfaceCoordinates('full_surface')

    vCoords = mesh.getCommonGrid()
    val = MPI.COMM_WORLD.reduce(numpy.sum(vCoords.flatten()),op=MPI.SUM)
    if MPI.COMM_WORLD.rank == 0:
        print('Sum of vCoords Inital:')
        reg_write(val,1e-8,1e-8)

    new_coords = coords0.copy()
    # Do a shearing sweep deflection:
    for i in xrange(len(coords0)):
        span = coords0[i,2]
        new_coords[i,0] += .05*span

    # Reset the newly computed surface coordiantes
    print('setting surface')
    mesh.setSurfaceCoordinates(new_coords, 'full_surface')
    print('warping mesh')
    mesh.warpMesh()
    
    # Get the sum of the warped coordinates
    #vCoords = mesh.getSolverGrid()
    vCoords = mesh.getWarpGrid()

    val = MPI.COMM_WORLD.reduce(numpy.sum(vCoords.flatten()),op=MPI.SUM)
    if MPI.COMM_WORLD.rank == 0:
        print('Sum of vCoords Warped:')
        reg_write(val,1e-8,1e-8)

    # Create a dXv vector to do test the mesh warping with:
    dXv_warp = numpy.linspace(0,1.0, mesh.warp.griddata.warpmeshdof)

    if not 'complex' in sys.argv:
        if MPI.COMM_WORLD.rank == 0:
            print('Computing Warp Deriv')
        mesh.warpDeriv(dXv_warp,solverVec=False)
        dXs = mesh.getdXs('full_surface')
            
        val = MPI.COMM_WORLD.reduce(numpy.sum(dXs.flatten()),op=MPI.SUM)
        if MPI.COMM_WORLD.rank == 0:
            print('Sum of dxs:')
            reg_write(val,1e-8,1e-8)
    else:
        
        # add a complex perturbation on all surface nodes simultaneously:
        for i in xrange(len(coords0)):
            new_coords[i,:] += h*1j

        # Reset the newly computed surface coordiantes
        mesh.setSurfaceCoordinates(new_coords, 'full_surface')
        mesh.warpMesh()  

        vCoords = mesh.getWarpGrid()
        deriv = numpy.imag(vCoords)/h
        deriv = numpy.dot(dXv_warp,deriv)
        val = MPI.COMM_WORLD.reduce(numpy.sum(deriv),op=MPI.SUM)
        val/=len(coords0)*3
        if MPI.COMM_WORLD.rank == 0:
            print('Sum of dxs:')
            reg_write(val,1e-8,1e-8)


    del mesh
    #os.system('rm -fr *.cgns *.dat')


if __name__ == '__main__':
    if len(sys.argv) == 1:
        test1()
        test2()
        test3()
        #test4()
        test5()
    else:
        # Run individual ones
        if 'test1' in sys.argv:
            test1()
        if 'test2' in sys.argv:
            test2()
        if 'test3' in sys.argv:
            test3()
        # if 'test4' in sys.argv:
        #     test4()
        if 'test5' in sys.argv:
            test5()
