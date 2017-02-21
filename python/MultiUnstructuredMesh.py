#!/usr/bin/python
"""
Multiple USMesh instances contained in this MultiUSMesh object.

The MultiUnstructuredMesh module is used for interacting with multiple
unstructured (or structured!) meshes - typically used in a 3D CFD
program.

It contains the following classes:

MultiUSMesh: General class for working with multiple unstructured meshes

Copyright (c) 2017 by John Jasa
All rights reserved. Not to be used for commercial purposes.
Revision: 1.0   $Date: 02/14/2017$

Developers:
-----------
- John Jasa (JPJ)
- Ney Secco (NRS)

History
-------
    v. 1.0 - Initial Class Creation (JPJ, 2017)
"""
from __future__ import print_function
from __future__ import division
# =============================================================================
# Imports
# =============================================================================
import os
import re
import shutil
import copy
import numpy as np
from pprint import pprint
from mpi4py import MPI
from .MExt import MExt

class Error(Exception):
    """
    Format the error message in a box to make it clear this
    was a expliclty raised exception.
    """
    def __init__(self, message):
        msg = '\n+'+'-'*78+'+'+'\n' + '| pyWarpUstruct Error: '
        i = 22
        for word in message.split():
            if len(word) + i + 1 > 78: # Finish line and start new one
                msg += ' '*(78-i)+'|\n| ' + word + ' '
                i = 1 + len(word)+1
            else:
                msg += word + ' '
                i += len(word)+1
        msg += ' '*(78-i) + '|\n' + '+'+'-'*78+'+'+'\n'
        print(msg)
        Exception.__init__(self)

class Warning(object):
    """
    Format a warning message
    """
    def __init__(self, message):
        msg = '\n+'+'-'*78+'+'+'\n' + '| pyWarpUstruct Warning: '
        i = 24
        for word in message.split():
            if len(word) + i + 1 > 78: # Finish line and start new one
                msg += ' '*(78-i)+'|\n| ' + word + ' '
                i = 1 + len(word)+1
            else:
                msg += word + ' '
                i += len(word)+1
        msg += ' '*(78-i) + '|\n' + '+'+'-'*78+'+'+'\n'
        print(msg)

# =============================================================================
# UnstructuredMesh class
# =============================================================================

class MultiUSMesh(object):
    """
    This is the main Unstructured Mesh. This mesh object is designed
    to interact with an structured or unstructured CFD solver though a
    variety of interface functions.
    """
    def __init__(self):
        """
        Create the MultiUSMesh object.
        """

        self.meshes = []

    def addMesh(self, mesh):
        self.meshes.append(mesh)

    def getSurfaceCoordinates(self):
        """Returns all defined surface coordinates on this processor

        Returns
        -------
        coords : numpy array size (N,3)
            Specified surface coordinates residing on this
            processor. This may be empty array, size (0,3)
        """
        #
        # self._setInternalSurface()
        # coords = np.zeros((self.nSurf, 3), self.dtype)
        # self.warp.getsurfacecoordinates(np.ravel(coords))

        for mesh in self.meshes:
            mesh.getSurfaceCoordinates()

        return coords

    def setSurfaceCoordinates(self, coordinates):
        """Sets all surface coordinates on this processor

        Parameters
        ----------
        coordinates : numpy array, size(N, 3)
            The coordinate to set. This MUST be exactly the same size as
            the array obtained from getSurfaceCoordinates()
        """
        # if len(coordinates) != self.nSurf:
        #     raise Error("Incorrect length of coordinates supplied to "
        #                 "setSurfaceCoordinates on proc %d. Expected "
        #                 "aray of length %d, received an array of length "
        #                 "%d."%(self.comm.rank, self.nSurf, len(coordinates)))
                # self.warp.setsurfacecoordinates(np.ravel(coordinates))

        for mesh in self.meshes:
            mesh.setSurfaceCoordinates(coordinates)


    def setExternalMeshIndices(self, ind):
        """
        Set the indices defining the transformation of an external
        solver grid to the original CGNS grid. This is required to use
        USMesh functions that involve the word "Solver" and warpDeriv.
        The indices must be zero-based.

        Parameters
        ----------
        ind : numpy integer array
            The list of indicies this processor needs from the common mesh file
        """

        for mesh in self.meshes:
            mesh.setExternalMeshIndices(ind)

        # self.warp.setexternalmeshindices(ind)

    def getSolverGrid(self):
        """Return the current grid in the order specified by
        setExternalMeshIndices(). This is the main routine for
        returning the deformed mesh to the external CFD solver.

        Returns
        -------
        solverGrid, numpy array, real: The resulting grid.
           The output is returned in flatted 1D coordinate
           format. The len of the array is 3*len(indices) as
           set by setExternalMeshIndices()

        """
        for mesh in self.meshes:
            solverGrid = mesh.getSolverGrid()
        # solverGrid = np.zeros(self.warp.griddata.solvermeshdof, self.dtype)
        # warpGrid = self.getWarpGrid()
        # self.warp.warp_to_solver_grid(warpGrid, solverGrid)

        return solverGrid

    def setSurfaceDefinition(self, pts, conn, faceSizes):
        """This is the master function that determines the definition of the
        surface to be used for the mesh movement. This surface may be
        supplied from an external solver (such as ADflow) or it may be
        generated by pyWarpUStruct internally.

        Parameters
        ----------
        pts : array, size (M, 3)
            Nodes on this processor
        conn : int array, size  (sum(faceSizes))
            Connectivity of the nodes on this processor
        faceSizes : int Array size (N)
            Treat the conn array as a flat list with the faceSizes giving
            connectivity offset for each element.
        """



        for mesh in self.meshes:
            mesh.setSurfaceDefinition(pts, conn, faceSizes)

    def getWarpGrid(self):
        """
        Return the current grid. This function is typically unused. See
        getSolverGrid for the more useful interface functionality.

        Returns
        -------
        warpGrid, numpy array, real: The resulting grid.
            The output is returned in flatted 1D coordinate
            format.
        """
        for mesh in self.meshes:
            volumeCoords = getWarpGrid()
        return volumeCoords
        # return self.warp.getvolumecoordinates(self.warp.griddata.warpmeshdof)

    def getdXs(self):
        """Return the current values in dXs. This is the result from a
        mesh-warp derivative computation.

        Returns
        -------
        dXs :  numpy array
            The specific components of dXs. size(N,3). This the same
            size as the array obtained with getSurfaceCoordiantes(). N may
            be zero if this processor does not have any surface coordinates.
        """
        # dXs = np.zeros((self.nSurf, 3), self.dtype)
        # self.warp.getdxs(np.ravel(dXs))
        #
        # return dXs

        for mesh in self.meshes:
            dXs = mesh.getdXs()

        return dXs

    def warpMesh(self):
        """
        Run the applicable mesh warping strategy.

        This will update the volume coordinates to match surface
        coordinates set with aceCoordinates()
        """

        for mesh in self.meshes:
            mesh.warpMesh()

    def warpDeriv(self, dXv, solverVec=True):
        """Compute the warping derivative (dXv/dXs^T)*Vec (where vec is the
        dXv argument to this function.

        This is the main routine to compute the mesh warping
        derivative.

        Parameters
        ----------
        solverdXv :  numpy array
            Vector of size external solver_grid. This is typically
            obtained from the external solver's dRdx^T * psi
            calculation.

        solverVec : logical
            Flag to indicate that the dXv vector is in the solver
            ordering and must be converted to the warp ordering
            first. This is the usual approach and thus defaults to
            True.

        Returns
        -------
        None. The resulting calculation is available from the getdXs()
        function.

        """
        if solverVec:
            dXvWarp = np.zeros(self.warp.griddata.warpmeshdof, self.dtype)
            self.warp.solver_to_warp_grid(dXv, dXvWarp)
        else:
            dXvWarp = dXv
        self.warp.warpderiv(dXvWarp)

    def warpDerivFwd(self, dXs, solverVec=True):
        """Compute the forward mode warping derivative

        This routine is not used for "regular" optimization; it is
        used for matrix-free type optimization. dXs is assumed to be
        the the peturbation on all the surface nodes.

        Parameters
        ----------
        dXs : array, size Nsx3
            This is the forward mode peturbation seed. Same size as the
            surface mesh from getSurfaceCoordinates().
        solverVec : bool
            Whether or not to convert to the solver ordering.

        Returns
        -------
        dXv : array
            The peturbation on the volume meshes. It may be
            in warp ordering or solver ordering depending on
            the solverVec flag.
        """

        dXvWarp = np.zeros(self.warp.griddata.warpmeshdof)
        self.warpMesh()
        self.warp.warpderivfwd(np.ravel(dXs), dXvWarp)
        if solverVec:
            dXv = np.zeros(self.warp.griddata.solvermeshdof)
            self.warp.warp_to_solver_grid(dXvWarp, dXv)
            return dXv
        else:
            return dXvWarp

    def verifyWarpDeriv(self, dXv=None, solverVec=True, dofStart=0,
                        dofEnd=10, h=1e-6, randomSeed=314):
        """Run an internal verification of the solid warping
        derivatives"""

        if dXv is None:
            np.random.seed(randomSeed) # 'Random' seed to ensure runs are same
            dXvWarp = np.random.random(self.warp.griddata.warpmeshdof)
        else:
            if solverVec:
                dXvWarp = np.zeros(self.warp.griddata.warpmeshdof, self.dtype)
                self.warp.solver_to_warp_grid(dXv, dXvWarp)
            else:
                dXvWarp = dXv

        self.warp.verifywarpderiv(dXvWarp, dofStart, dofEnd, h)

# ==========================================================================
#                        Output Functionality
# ==========================================================================
    def writeGrid(self, fileName=None):
        """
        Write the current grid to the correct format

        Parameters
        ----------
        fileName : str or None
            Filename for grid. Should end in .cgns for CGNS files. It is
            not required for openFOAM meshes. This call will update the
            'points' file.
        """

        if self.meshType.lower() == 'cgns':
            # Copy the default and then write
            if self.comm.rank == 0:
                shutil.copy(self.solverOptions['gridFile'], fileName)
            self.warp.writecgns(fileName)

        elif self.meshType.lower() == 'openfoam':
            self._writeOpenFOAMVolumePoints(self.getCommonGrid())

# =========================================================================
#                     Internal Private Functions
# =========================================================================

    def _setSymmetryConditions(self):
        """This function determines the symmetry planes used for the
        computation. It has a similar structure to setInternalSurface.

        Symmetry planes are specified throught 'symmetryPlanes' option, which
        has the form

        'symmetryPlanes':[[[x1,y1, z1],[dir_x1, dir_y1, dir_z1]],[[x2,y2, z2],[dir_x2, dir_y2, dir_z2]],...]

        Examples
        --------
        meshOptions = {'symmetryPlanes':[[[0.,0., 0.],[0., 1., 0.]]]}
        mesh = USMesh(options=meshOptions,comm=gcomm)

        """

        symmList = self.solverOptions['symmetryPlanes']
        if symmList is not None:
            # The user has explictly supplied symmetry planes. Use those
            pts = []
            normals = []
            for i in range(len(symmList)):
                pts.append(symmList[i][0])
                normals.append(symmList[i][1])

        else:
            # Otherwise generate from the geometry.
            planes = []
            if self.meshType.lower() == "cgns":
                if self.comm.rank == 0:

                    # Do the necessary fortran preprocessing
                    if self.warp.cgnsgrid.cgnsstructured:
                        self.warp.processstructuredpatches()
                    else:
                        self.warp.processunstructuredpatches()

                    fullConn = self.warp.cgnsgrid.surfaceconn-1
                    fullPts = self.warp.cgnsgrid.surfacepoints
                    fullPatchNames = self._processFortranStringArray(
                        self.warp.cgnsgrid.surfacenames)

                    symmFamilies = set()
                    if self.solverOptions['symmetrySurfaces'] is None:
                        # Use all symmetry surfaces:
                        for i in range(len(self.warp.cgnsgrid.surfaceissymm)):
                            if self.warp.cgnsgrid.surfaceissymm[i]:
                                symmFamilies.add(fullPatchNames[i].lower())
                    else:
                        # The user has supplied a list of surface families
                        for name in self.solverOptions['symmetrySurfaces']:
                            symmFamilies.add(name.lower())

                usedFams = set()
                if self.warp.cgnsgrid.cgnsstructured:
                    if self.comm.rank == 0:

                        # Pull out data and convert to 0-based ordering
                        fullPatchSizes = self.warp.cgnsgrid.surfacesizes.T

                        # Now we loop back over the "full" versions of
                        # thing and just take the ones that correspond
                        # to the families we are using.

                        curNodeIndex = 0
                        curCellIndex = 0
                        curOffset = 0
                        for i in range(len(fullPatchNames)):
                            curNodeSize = fullPatchSizes[i][0]*fullPatchSizes[i][1]
                            curCellSize = (fullPatchSizes[i][0]-1)*(fullPatchSizes[i][1]-1)

                            if fullPatchNames[i] in symmFamilies:
                                # Keep track of the families we've actually used
                                usedFams.add(fullPatchNames[i])

                                # Determine the average normal for these points:
                                conn =  fullConn[curCellIndex:curCellIndex+curCellSize*4] - fullConn[curCellIndex] + 1
                                pts =   fullPts[curNodeIndex:curNodeIndex+curNodeSize*3]
                                avgNorm = self.warp.averagenormal(
                                    pts, conn, 4*np.ones(curCellSize, 'intc'))
                                planes.append([pts[0:3], avgNorm])
                            else:
                                # If we skipped, we increment the offset
                                curOffset += curNodeSize

                            curNodeIndex += curNodeSize*3
                            curCellIndex += curCellSize*4

                    # end for (root proc)

                else: # unstructured

                    # We won't do this in general. The issue is that
                    # each of the elements needs to be checked
                    # individually since one sym BC may have multiple
                    # actual symmetry planes. This could be done, but
                    # since plane elemination code is in python and
                    # slow, we'll just defer this and make the user
                    # supply the symmetry planes.

                    raise Error("Automatic determine of symmetry surfaces is "
                                "not supported for unstructured CGNS. Please "
                                "specify the symmetry planes using the "
                                "'symmetryPlanes' option. See the _setSymmetryConditions()"
                               " documentation string for the option prototype.")

                # Check if all supplied family names were actually
                # used. The user probably wants to know if a family
                # name was specified incorrectly.
                if self.comm.rank == 0:
                    if usedFams < symmFamilies:
                        missing = list(symmFamilies - usedFams)
                        Warning("Not all specified symm families that "
                                "were given were found the CGNS file. "
                                "The families not found are %s."%(repr(missing)))

            elif self.meshType.lower() == "openfoam":

                # We could probably implement this at some point, but
                # it is not critical

                raise Error("Automatic determine of symmetry surfaces is "
                            "not supported for OpenFoam meshes. Please "
                            "specify the symmetry planes using the "
                            "'symmetryPlanes' option. See the _setSymmetryConditions()"
                            " documentation string for the option prototype.")

            # Now we have a list of planes. We have to reduce them to the
            # set of independent planes. This is tricky since you can have
            # have to different normals belonging to the same physical
            # plane. Since we don't have that many, we just use a dumb
            # double loop.

            def checkPlane(p1, n1, p2, n2):
                # Determine if two planes defined by (pt, normal) are
                # actually the same up to a normal sign.

                # First check the normal...if these are not the same,
                # cannot be the same plane
                if abs(np.dot(n1, n2)) < 0.99:
                    return False

                # Normals are the same direction. Check if p2 is on the
                # first plane up to a tolerance.

                d = p2 - p1
                d1 = p2 - np.dot(d, n1)*n1

                if np.linalg.norm(d1 - p2) / (np.linalg.norm(d) + 1e-30) > 1e-8:
                    return False

                return True

            uniquePlanes = []
            flagged = np.zeros(len(planes), 'intc')
            for i in range(len(planes)):
                if not flagged[i]:
                    uniquePlanes.append(planes[i])
                    curPlane = planes[i]
                    # Loop over remainder to check:
                    for j in range(i+1, len(planes)):
                        if checkPlane(curPlane[0], curPlane[1],
                                      planes[j][0], planes[j][1]):
                            flagged[j] = 1

            # Before we return, reset the point for each plane to be as
            # close to the origin as possible. This will help slightly
            # with the numerics.

            pts = []
            normals = []

            for i in range(len(uniquePlanes)):
                p = uniquePlanes[i][0]
                n = uniquePlanes[i][1]

                p2 = np.zeros(3)
                d = p2 - p
                pstar = p2 - np.dot(d, n)*n

                normals.append(n)
                pts.append(pstar)

        normals = self.comm.bcast(np.array(normals))
        pts = self.comm.bcast(np.array(pts))

        # Let the user know what they are:
        if self.comm.rank == 0:
            print ('+-------------------- Symmetry Planes -------------------+')
            print ('|           Point                        Normal          |')
            for i in range(len(pts)):
                print ('| (%7.3f %7.3f %7.3f)    (%7.3f %7.3f %7.3f) |'%(
                    np.real(pts[i][0]), np.real(pts[i][1]), np.real(pts[i][2]),
                    np.real(normals[i][0]), np.real(normals[i][1]), np.real(normals[i][2])))
            print ('+--------------------------------------------------------+')

        # Now set the data into fortran.
        self.warp.setsymmetryplanes(pts.T, normals.T)

# =========================================================================
#                     Utility functions
# =========================================================================

    def _setMeshOptions(self):
        """ Private function to set the options currently in
        self.solverOptions to the corresponding place in Fortran"""
        o = self.solverOptions
        self.warp.gridinput.ldeffact = o['LdefFact']
        self.warp.gridinput.alpha = o['alpha']
        self.warp.gridinput.aexp = o['aExp']
        self.warp.gridinput.bexp = o['bExp']
        self.warp.gridinput.symmtol = o['symmTol']
        self.warp.gridinput.userotations = o['useRotations']
        self.warp.gridinput.zerocornerrotations = o['zeroCornerRotations']
        self.warp.gridinput.cornerangle = o['cornerAngle']
        self.warp.gridinput.errtol = o['errTol']
        self.warp.kd_tree.bucket_size = o['bucketSize']
        if o['evalMode'].lower() == 'fast':
            self.warp.gridinput.evalmode = self.warp.gridinput.eval_fast
        elif o['evalMode'].lower() == 'exact':
            self.warp.gridinput.evalmode = self.warp.gridinput.eval_exact

    def _checkOptions(self, solverOptions):
        """ Check the solver options against the default ones and add
        option iff it is NOT in solverOptions """

        for key in self.solverOptionsDefault.keys():
            if not key in solverOptions.keys():
                solverOptions[key] = self.solverOptionsDefault[key]
            else:
                self.solverOptionsDefault[key] = solverOptions[key]

        # Print a couple of warnings about aExp and bExp which are not
        # fully implemented.
        if self.comm.rank == 0:
            if self.solverOptions['aExp'] != 3.0 or self.solverOptions['bExp'] != 5.0:
                Warning('The aExp and bExp options are currently not implemented '
                        'and will always use their default values of aExp=3.0 and '
                        'bExp=5.0')

        return solverOptions

    def _printCurrentOptions(self):
        """Prints a nicely formatted dictionary of all the current options to
        the stdout on the root processor
        """
        if self.comm.rank == 0:
            print('+---------------------------------------+')
            print('|     All pyWarpUstruct Options:        |')
            print('+---------------------------------------+')
            pprint(self.solverOptions)

    def _processFortranStringArray(self, strArray):
        """Getting arrays of strings out of Fortran can be kinda nasty. This
        takes the array and returns a nice python list of strings"""
        shp = strArray.shape
        arr = strArray.reshape((shp[1],shp[0]), order='F')
        tmp = []
        for i in range(arr.shape[1]):
            tmp.append(''.join(arr[:, i]).strip().lower())

        return tmp

    def __del__(self):
        """Release all the mesh warping memory. This should be called
        automatically when the object is garbage collected."""
        for mesh in self.meshes:
            mesh.warp.releasememory()
