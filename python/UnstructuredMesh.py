#!/usr/bin/python
"""Unstructured Mesh

The UnstructuredMesh module is used for interacting with an
unstructured (or structured!)  mesh - typically used in a 3D CFD
program.

It contains the following classes:

USMesh: General class for working with unstructured meshes

Copyright (c) 2014 by C.A.(Sandy) Mader
All rights reserved. Not to be used for commercial purposes.
Revision: 1.0   $Date: 09/07/2014$

Developers:
-----------
- C.A.(Sandy) Mader (CAM)
- Gaetan. K. W. Kenway (GKK)
History
-------
	v. 1.0 - Initial Class Creation (CAM, 2014)
"""
from __future__ import print_function
from __future__ import division
# =============================================================================
# Imports
# =============================================================================
import os
import re
import shutil
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

# =============================================================================
# UnstructuredMesh class
# =============================================================================
IWALL = 0
IFAR = 1

class USMesh(object):
    """
    This is the main Unstructured Mesh. This mesh object is designed to
    interact with an unstructured CFD solver though a variety of
    interface functions.
 
    """
    def __init__(self, options=None, comm=None, debug=False):
        """
        Create the USMesh object. 

        Parameters
        ----------
        options :  dictionary
            A dictionary containing the options for the mesh movement 
            strategies. THIS IS NOT OPTIONAL.

        comm : MPI_INTRA_COMM
            MPI communication (as obtained from mpi4py) on which to create 
            the USMesh object. If not provided, MPI_COMM_WORLD is used. 

        debug : bool
            Flag specifying if the MExt import is automatically deleted. 
            This needs to be true ONLY when a symbolic debugger is used.
        """

        if comm is None:
            self.comm = MPI.COMM_WORLD
        else:
            self.comm = comm

        # Check if warp has already been set if this is this has been
        # interhited to complex version
        try: 
            self.warp
        except AttributeError:
            curDir = os.path.dirname(os.path.realpath(__file__))
            self.warp = MExt('warpustruct', [curDir], debug=debug)._module
        if options is not None:
            self.solverOptions = options
        else:
            raise Error("The 'options' keyword argument is *NOT* "
                        "optional. An options dictionary must be passed upon "
                        "creation of this object")

        # Initialize PETSc if not done so
        self.warp.initpetsc(self.comm.py2f())

        # Set realtype of 'd'. 'D' is used in Complex and set in
        # UnstructuredMesh_C.py
        self.dtype = 'd'

        # Defalut options for mesh warping
        self.solverOptionsDefault = {
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
        }
        
        # Set remaining default values
        self._checkOptions(self.solverOptions)
        self._setMeshOptions()
        self.printCurrentOptions()

        # Initialize various bits of stored information
        self.OFData = {}
        self.warpInitialized = False
        self.patchSizes = None
        self.familyList = None
        self.patchPtr = None
        self.patchNames = None
        self.faceSizes = None
        self.conn = None
        fileName = self.solverOptions['gridFile']

        # Determine how to read:
        if os.path.isfile(fileName):
            self.warp.readstructuredcgns(fileName)
            self._getCGNSFamilyList()

        elif os.path.isdir(fileName):
            self._readOFGrid(fileName)
        else:
            raise Error('Invalid input file type. CGNS or OpenFOAM required')

        # Add a default family group 'all' which will always be there.
        self.familyGroup = {}
        self.addFamilyGroup('all')

# =========================================================================
#                  Local mesh reading functionality
# =========================================================================

    def _readOFGrid(self, caseDir):
        """Read in the mesh points and connectivity from the polyMesh
        directory
    
        Parameters
        ----------
        caseDir : str
            The directory containing the openFOAM Mesh files
        """
        
        # Copy the reference points file to points to ensure
        # consistant starting point
        self._getOFFileNames(caseDir)
        shutil.copyfile(self.OFData['refPointsFile'], self.OFData['pointsFile'])
      
        # Read in the volume points
        self._readOFVolumeMeshPoints()

        # Read the face info for the mesh
        self._readOFFaceInfo()

        # Read the boundary info
        self._readOFBoundaryInfo()
      
        # Read the cell info for the mesh
        self._readOFCellInfo()

        # Finally create the underlying data structure:
        faceSizes = []
        faceConn = []
        pts = []
        x0 = self.OFData['x0']
        faces = self.OFData['faces']
        connCount = 0
        localFamilies = []
        localType = []
        self.patchNames = []
        self.patchPtr = [0]
        symNodes = []
        for bName in self.OFData['boundaries']:
            bType = self.OFData['boundaries'][bName]['type'].lower()
            if bType in ['patch', 'wall']:
                localFamilies.append(bName)
                if bType == 'wall':
                    localType.append(IWALL)
                if bType == 'patch':
                    localType.append(IFAR)

                # Apparently openfoam will list boundaries with zero
                # faces on them:
                nFace = len(self.OFData['boundaries'][bName]['faces'])
                if nFace > 0:
                    for iFace in self.OFData['boundaries'][bName]['faces']:
                        face = faces[iFace]
                        nNodes = len(face)
                        # Pull out the 'len(face)' nodes from x0.
                        for j in range(nNodes):
                            pts.append(x0[face[j]])
                        faceSizes.append(nNodes)
                        faceConn.extend(range(connCount, connCount + nNodes))
                        connCount += nNodes

                    self.patchNames.append(bName.lower())
                    self.patchPtr.append(connCount)
            if bType.lower() == 'symmetry':
                # Gather up all the sym nodes
                for iFace in self.OFData['boundaries'][bName]['faces']:
                    face = faces[iFace]
                    for j in range(len(face)):
                        symNodes.append(x0[face[j]])

        # Gather all symmetry nodes:
        if len(symNodes) == 0:
            symNodes = np.zeros((0, 3))
        else:
            symNodes = np.array(symNodes)
        
        allSymNodes = self.comm.allgather(symNodes)
        allSymNodes = np.vstack(allSymNodes)
        iSymm = 0
        if len(allSymNodes > 0):
            symmSum = np.sum(allSymNodes, axis=0)
            n = len(symmSum)
            if symmSum[0] / n < self.solverOptions['symmTol']:
                iSymm = 1
            elif symmSum[1] / n < self.solverOptions['symmTol']:
                iSymm = 2
            elif symmSum[2] / n < self.solverOptions['symmTol']:
                iSymm = 3
            self.warp.gridinput.isymm = iSymm

        # Create the internal structures for volume mesh
        wallNodes = np.zeros(len(x0), 'intc')
        self.warp.createcommongrid(x0.T, wallNodes)
        pts = np.array(pts)
        # And run the common routine for setting up the surface
        self.warp.initializewarping(np.ravel(pts.real.astype('d')),
                                    faceSizes, faceConn)
        self.warpInitialized = True

        # Now communicate all the family stuff
        tmpFamilies = self.comm.allgather(localFamilies)
        tmpPatchType = self.comm.allgather(localType)
        
        self.familyList = {}
        for iProc in range(self.comm.size):
            for j in range(len(tmpFamilies[iProc])):
                fam = tmpFamilies[iProc][j]
                pType = tmpPatchType[iProc][j]
                if fam not in self.familyList:
                    self.familyList[fam.lower()] = pType

    def _readOFVolumeMeshPoints(self):
        '''
        Return an numpy array of the mesh points
        '''

        # Open the points file for reading
        f = open(self.OFData['pointsFile'], 'r')
        lines = f.readlines()
        f.close()
        
        # Read Regular Header
        foamHeader, i, N = self._parseFoamHeader(lines, i=0)
        
        try:
            fmt = foamHeader['format'].lower() 
        except:
            fmt = 'ascii' # assume ascii if format is missing.

        x = np.zeros((N, 3), self.dtype)
        pointLine = re.compile(r'\((\S*)\s*(\S*)\s*(\S*)\)\s*\n')

        if fmt == 'binary':
            raise Error('Binary Reading is not yet supported.')
        else:
            k = 0
            for j in range(N):
                res = pointLine.match(lines[j + i])
                if res:
                    k += 1
                    for idim in range(3):
                        x[j, idim] = float(
                            lines[j + i][res.start(idim+1):res.end(idim+1)])

            if k != N:
                raise Error('Error reading grid coordinates. Expected %d'
                            ' coordinates but found %d'%(N, k))
        self.OFData['x0'] = x

    def _readOFBoundaryInfo(self):
        '''
        read the boundary file information for this case and store in a dict.
        '''

        # Open the boundary file
        f = open(self.OFData['boundaryFile'], 'r')
        lines = f.readlines()
        f.close()

        # Read Regular Header
        foamHeader, i, N = self._parseFoamHeader(lines, i=0)

        # We don't actually need to know the number of blocks...just
        # read them all:
        boundaries = {}
        keyword = re.compile(r'\s*([a-zA-Z]{1,100})\s*\n')
        while i < len(lines):
            res = keyword.match(lines[i])
            if res:
                boundaryName = lines[i][res.start(1):res.end(1)]
                i, blkData = self._readBlock(lines, i)

                # Now we can just read out the info we need:
                startFace = int(blkData['startFace'])
                nFaces = int(blkData['nFaces'])
                faces = np.arange(startFace, startFace+nFaces)
                nodes = []
                for face in faces:
                    nodes.extend(self.OFData['faces'][face])
                boundaries[boundaryName] = {
                    'faces':faces,
                    'nodes':nodes,
                    'type':blkData['type']}
            else:
                i += 1
        self.OFData['boundaries'] = boundaries

    def _readOFFaceInfo(self):
        """
        Read the face info for this case.
        """

        f = open(self.OFData['faceFile'], 'r')
        lines = f.readlines()
        f.close()
        
        # Read Regular Header
        foamHeader, i, N = self._parseFoamHeader(lines, i=0)
        
        try:
            fmt = foamHeader['format'].lower() 
        except:
            fmt = 'ascii' # assume ascii if format is missing.

        if fmt == 'binary':
            raise Error('Binary Reading is not yet supported.')
        else:
            faces = [[] for k in range(N)]
            for j in range(N):
                line = lines[j+i].replace('(', ' ')
                line = line.replace(')', ' ')
                aux = line.split()
                nNode = int(aux[0])
                for k in range(1, nNode+1):
                    faces[j].append(int(aux[k]))
        
            self.OFData['faces'] = faces

    def _readOFCellInfo(self):
        """Read the boundary file information for this case and store in the
        OFData dictionary."""

        # ------------- Read the owners file --------------- 
        f = open(self.OFData['ownerFile'], 'r')
        lines = f.readlines()
        f.close()

        # Read Regular Header
        foamHeader, i, N = self._parseFoamHeader(lines, i=0)
        
        try:
            fmt = foamHeader['format'].lower() 
        except:
            fmt = 'ascii' # assume ascii if format is missing.

        owner = np.zeros(N, 'intc')

        if fmt == 'binary':
            raise Error('Binary Reading is not yet supported.')
        else:
            for j in range(N):
                owner[j] = int(lines[j + i])

        # ------------- Read the neighbour file --------------- 
        f = open(self.OFData['neighbourFile'], 'r')
        lines = f.readlines()
        f.close()

        # Read Regular Header
        foamHeader, i, N = self._parseFoamHeader(lines, i=0)
        
        try:
            fmt = foamHeader['format'].lower() 
        except:
            fmt = 'ascii' # assume ascii if format is missing.

        neighbour = np.zeros(N, 'intc')
        if fmt == 'binary':
            raise Error('Binary Reading is not yet supported.')
        else:
            for j in range(N):
                neighbour[j] = int(lines[j + i])

        self.OFData['owner'] = owner
        self.OFData['neighbour'] = neighbour

    def addFamilyGroup(self, groupName, families=None):
        """Create a grouping of CGNS families called "groupName"
        
        Parameters
        ----------
        groupName : str
            The name to call this collection of families
        families : list
           An (optional) list containing the names of the CGNS familes
           the user wants included in "groupName". The items in the
           list must be valid family names. The user can call
           printFamilyList() to determine what families are present in
           the CGNS grid.  Default: All families
        """
      
        # Use whole list by default
        if families == None: 
            families = []
            for fam in self.familyList:
                if self.familyList[fam] == IWALL:
                    families.append(fam)
                    
        groupFamList = []
        for fam in families:
            if fam.lower() in self.familyList:
                if self.familyList[fam.lower()] == IWALL:
                    groupFamList.append(fam.lower())
            else:
                raise Error("%s was not in family list!"% fam)
        
        self.familyGroup[groupName] = {'families':groupFamList}

        if self.warpInitialized: # Warp is initialized so we can continue
            nodeIndices = []
            nodeCount = 0

            # Loop over all possible patch Names
            for i in range(len(self.patchNames)):
                sizeOfPatch = self.patchPtr[i+1] - self.patchPtr[i]
                
                # Check if the current patch is part of what we are adding
                if self.patchNames[i] in groupFamList:
                    # Get the nodes indices:
                    nodeIndices.extend(range(
                        3*self.patchPtr[i], 3*self.patchPtr[i+1]))
                nodeCount += sizeOfPatch

            self.familyGroup[groupName]['indices'] = nodeIndices

            if self.patchSizes is not None:
                # Now we will compute the connectivity of the
                # requested family group. This is currently only
                # possible when the mesh object has been initialized
                # using and external structured surface mesh.
                cellConn = []
                nNodesSkipped = 0
                cCnt = 0
          
                for i in range(len(self.patchNames)):
                    cellSize = (self.patchSizes[i][0]-1)*\
                               (self.patchSizes[i][1]-1)
                    if self.patchNames[i] in groupFamList: 
                        # Add the connecitivity minus the number of nodes
                        # we've skipped so far. 
                        cellConn.extend(
                            self.conn[4*cCnt:4*(cCnt+cellSize)].copy() - 
                            nNodesSkipped)
                    else:
                        # If we didn't add this patch increment the number
                        # of nodes we've skipped
                        nNodesSkipped += (self.patchPtr[i+1] - self.patchPtr[i])
                    cCnt += cellSize

                self.familyGroup[groupName]['connectivity'] = (
                    np.array(cellConn).astype('intc').reshape((len(cellConn)/4, 4)))

    def printFamilyList(self):
        """Prints the families in the CGNS file"""
        if self.comm.rank == 0:
            print('Family list is:', self.familyList.keys())

    def getSurfaceCoordinates(self, groupName):
        """ 
        Returns all surface coordinates on this processor in group
        'groupName'

        Parameters
        ----------
        groupName : str
           The group from which to obtain the coordinates.
           This name must have been obtained from addFamilyGroup() or
           be the default 'all' which contains all surface coordinates

        Returns
        -------
        coords : numpy array size (N,3)
            Coordinates of the requested group. This may be empty 
            array, size (0,3)
        """
        self._setInternalSurface()
        indices = self._getIndices(groupName)
        coords = np.zeros((len(indices)//3, 3), self.dtype)
        self.warp.getsurfacecoordinates(indices, np.ravel(coords))

        return coords
   
    def getSurfaceConnectivity(self, groupName):
        """
        Returns the connectivities of the surface coordinates that were
        obtained with getSurfaceCoordinates()
        """
        return self.familyGroup[groupName]['connectivity'].flatten()

    def setSurfaceCoordinates(self, coordinates, groupName):
        """Sets all surface coordinates on this processor in group
        "groupName"

        Parameters
        ----------
        coordinates : numpy array, size(N, 3)
            The coordinate to set. This MUST be exactly the same size as 
            the array obtained from getSurfaceCoordinates()

        groupName : str
            The group to set the coordinate in. This name must have
            been set using addFamilyGroup() or be the default
            'all' which contains all surface coordiantes 
        """

        indices = self._getIndices(groupName)
        self.warp.setsurfacecoordinates(indices, np.ravel(coordinates))

    def setExternalMeshIndices(self, ind):
        """ 
        Set the indicies defining the transformation of an external
        solver grid to the original CGNS grid. This is required to use
        USMesh functions that involve the word "Solver" and warpDeriv
        
        Parameters
        ----------
        ind : numpy integer array
            The list of indicies this processor needs from the common mesh file
        """
        self.warp.setexternalmeshindices(ind)

    def setExternalSurface(self, patchNames, patchSizes, conn, pts):
        """
        This is the Master function that defines the surface
        coordinates the external solver wants to use.

        This function is **not** typically called by the user, but is
        done automatically by an external solver (ex. SUmb) in the
        setMesh() call. 

        Parameters
        ----------
        patchName : list of strings, length nPatch
            The names of each patch
        patchSizes : array size (nPatch, 2)
            The size of sizes of the patches
        conn : int array, size (N, 4)
            Connectivity of the nodes on this processor
        pts : array, size (M, 3)
            Nodes on this processor. 
        """

        # Since this routine is called by an external structured
        # solver we will unstructured-itize the data. Essentially we
        # just turn the patch sizes in to pointer of length patchNames
        # + 1 such that we can figure out what indices correspond to
        # the particular family. 

        patchPtr = np.zeros(len(patchSizes) + 1, 'intc')
        for i in range(len(patchSizes)):
            patchPtr[i+1] = patchPtr[i] + patchSizes[i][0]*patchSizes[i][1]
        self.patchSizes = patchSizes
        self.patchPtr = patchPtr
        self.patchNames = patchNames
        self.conn = np.array(conn).flatten()
        # All faces are quads so just use '4' for them all
        self.faceSizes = 4*np.ones(len(conn), 'intc') 
        
        # Call the fortran initialze warping command with  the
        # coordinates and the patch connectivity given to us. 
        self.warp.initializewarping(np.ravel(pts.real.astype('d')),
                                    self.faceSizes, self.conn)
        self.warpInitialized = True
        
        # We can now back out the indices that should go along with
        # the groupNames that may have already been added. 
        for key in self.familyGroup.keys():
            self.addFamilyGroup(key, self.familyGroup[key]['families'])

    def getSolverGrid(self):
        """
        Return the current grid in the order specified by
        setExternalMeshIndices()

        Returns
        -------
        solverGrid, numpy array, real: The resulting grid. 
           The output is returned in flatted 1D coordinate
           format. The len of the array is 3*len(indices) as
           set by setExternalMeshIndices()
        """
        solverGrid = np.zeros(self.warp.griddata.solvermeshdof, self.dtype)
        warpGrid = self.getWarpGrid()
        self.warp.warp_to_solver_grid(warpGrid, solverGrid)
        
        return solverGrid

    def getWarpGrid(self):
        """
        Return the current grid. This funtion is typically unused. See
        getSolverGrid for the more useful interface functionality. 
        
        Returns
        -------
        warpGrid, numpy array, real: The resulting grid. 
            The output is returned in flatted 1D coordinate
            format. 
        """
        return self.warp.getvolumecoordinates(self.warp.griddata.warpmeshdof)
        
    def getCommonGrid(self):
        """Return the grid int he original ordering. This is required for the
        openFOAM tecplot writer since the connectivity is only known
        in this ordering. 
        """
        return self.warp.getcommonvolumecoordinates(
            self.warp.griddata.commonmeshdof)

    def getdXs(self, groupName):

        """
        Return the current values in dXs. This is typically the result
        of a mesh warp-derivative computation.

        Parameters
        ----------
        groupName : str
            The group from which to obtain the derivatives.  This name
            must have been obtained from addFamilyGroup() or be the
            default 'all' which contains all surface coordiantes

        Returns
        -------
        dXs :  numpy array
            The specific components of dXs.  size(N,3) where N is the
            number of points associated with the families defined by
            'groupName'
            """

        indices = self._getIndices(groupName)
        dXs = np.zeros((len(indices)//3, 3), self.dtype)
        self.warp.getdxs(indices, np.ravel(dXs))

        return dXs

    def setdXs(self, dXs, groupName):
        """
        Set just the components into dXs given by groupName

        Parameters
        ----------
        dXs : numpy array
            The components correponding to groupName to set in dXs.

        groupName : str
            The group defining where components will be set in dXs
            """
        
        indices = self._getIndices(groupName)
        self.warp.setdxs(indices, np.ravel(dXs))

    def sectionVectorByFamily(self, groupName, inVec):
        """
        Partition a dXs-like vector (inVec) that is 'full size', ie
        size of the 'all' group and return just the part corresponding
        to 'groupName'.

        Parameters
        ----------
        groupName : str
            The group defining where components will be returned
        inVec : real numpy array
            The components to section whose length corresponds to the
            len of the 'all' group groupName to set in dXs.

        Returns
        -------
        dXs : numpy array
            The section of inVec that corresponds to the groupName
            """
        indices = self._getIndices('all')
        self.warp.setdxs(indices, np.ravel(inVec))
        indices = self._getIndices(groupName)
        outVec = np.zeros((len(indices)//3, 3), self.dtype)
        self.warp.getdxs(indices, np.ravel(outVec))

        return outVec

    def expandVectorByFamily(self, groupName, inVec):
        """
        Expand a dXs like vector (inVec) that correponds to
        'groupName' and return a vector that is of full size, ie the
        'all' group. It does the reverse operation of
        :func:`sectionVectorByFamily`. All other entries of the output
        are zeroed.

        Parameters
        ----------
        groupName : str
            The group defining where components will be added
        inVec : numpy array
            The components to section whose length corresponds to the
            len of 'groupName' section

        Returns
        -------
        dXs : numpy array
            The expanded form of dXs corresponding to the 'all' group
            """
        indices = self._getIndices(groupName)
        self.warp.setdxs(indices, np.ravel(inVec))
        indices = self._getIndices('all')
        outVec = np.zeros((len(indices)//3, 3), self.dtype)
        self.warp.getdxs(indices, np.ravel(outVec))
        
        return outVec

    def warpMesh(self):
        """ 
        Run the applicable mesh warping strategy.

        This will update the volume coordinates to match surface
        coordinates set with setSurfaceCoordinates()
        """
        self._setMeshOptions()
        self.warp.warpmesh()

    def warpDeriv(self, dXv, solverVec=True, surfOnly=False):
        """Compute the warping derivative (dXvdXs^T)*Vec

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

        surfOnly : logical
            Flag to indicate that a "fake" mesh warp is to be done. In
            this case, only dXv values that are ALREADY on the surface
            are taken. This is only used in a few select cases and
            thus defaults to False.

        Returns
        -------
        None. However, the resulting calculation is available from
        the getdXs() function.
        """

        if not surfOnly:
            if solverVec:
                dXvWarp = np.zeros(self.warp.griddata.warpmeshdof, self.dtype)
                self.warp.solver_to_warp_grid(dXv, dXvWarp)
            else:
                dXvWarp = dXv
            self.warp.warpderiv(dXvWarp)
        else:
            if not solverVec:
                raise Error('Fake warpDeriv only works with solverVec.')
            self.warp.warpderivsurfonly(dXv)

    def verifyWarpDeriv(self, dXv=None, solverVec=True, dofStart=0, 
                        dofEnd=10, h=1e-6):
        """Run an internal verification of the solid warping
        derivatives"""
        
        if dXv is None:
            np.random.seed(314) # 'Random' seed to ensure runs are same
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
    def writeGridCGNS(self, fileName):
        """
        Write the current grid to a CGNS file
        """
        raise Error("Not implemented yet")
        
    def writeOFGridTecplot(self, fileName):
        """
        Write the current openFOAM grid to a tepplot FE polyhedron file

        Parameters
        ----------
        fileName : str
            Filename to use. Should end in .dat for tecplot ascii file.
        """
        if not self.OFData:
            return
        if self.comm.size == 1:
            f = open(fileName, 'w')
        else:
            fname, fext = os.path.splitext(fileName)
            fileName = fname + '%d'%self.comm.rank + fext
            f = open(fileName, 'w')

        # Extract the data we need from the OFDict to make the code a
        # little easier to read:
        nodes = self.getCommonGrid()
        nodes = nodes.reshape((len(nodes)/3, 3))
        nPoints = len(nodes)

        faces = self.OFData['faces']
        nFaces = len(faces)

        owner = self.OFData['owner']
        nCells = np.max(owner) + 1
        
        neighbour = self.OFData['neighbour']
        nNeighbours = len(neighbour)

        # write the node numbers for each face
        faceNodeSum = 0
        for i in range(nFaces):   
            faceNodeSum += len(faces[i])

        # Tecplot Header Info:
        f.write('TITLE = "Example Grid File"\n')
        f.write('FILETYPE = GRID\n')
        f.write('VARIABLES = "X" "Y" "Z"\n')
        f.write('Zone\n')
        f.write('ZoneType=FEPOLYHEDRON\n')
        f.write('NODES=%d\n'%nPoints)
        f.write('FACES=%d\n'%nFaces)
        f.write('ELEMENTS=%d\n'%nCells)
        f.write('TotalNumFaceNodes=%d\n'%faceNodeSum)
        f.write('NumConnectedBoundaryFaces=%d\n'%0)
        f.write('TotalNumBoundaryConnections=%d\n'%0)

        # Write the points to file in block data format
        for idim in range(3):
            for j in range(nPoints):
                f.write('%g\n'% nodes[j, idim])

        # Write the number of nodes in each face
        for i in range(nFaces):
            f.write('%d '%len(faces[i]))
            if i%300 == 0:
                f.write('\n')

        f.write('\n')
        # write the node numbers for each face (+1 for tecplot 1-based
        # ordering)
        for i in range(nFaces):   
            nPointsFace = len(faces[i])
            for j in range(nPointsFace):
                f.write('%d '% (faces[i][j]+1))
            f.write('\n')

        # write left elements (+1 for 1 based ordering)
        for i in range(nFaces):
            f.write('%d\n'%(owner[i]+1))

        # write right elements (+1 for 1 based ordering)
        for i in range(nFaces):
            if i < nNeighbours:
                f.write('%d\n'% (neighbour[i]+1))
            else:
                f.write('%d\n'%(0))
                
        f.close()

    def writeOpenFOAMVolumePoints(self):
        '''
        Write the most recent points to a file
        '''
        fileName = self.OFData['pointsFile']
        f = open(fileName, 'w')

        # write the file header
        f.write('/*--------------------------------*- C++ -*----------------------------------*\ \n')
        f.write('| =========                 |                                                 |\n')
        f.write('| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n')
        f.write('|  \\\\    /   O peration     | Version:  2.3.x                                 |\n')
        f.write('|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |\n')
        f.write('|    \\\\/     M anipulation  |                                                 |\n')
        f.write('\*---------------------------------------------------------------------------*/\n')
        f.write('FoamFile\n')
        f.write('{\n')
        f.write('    version     2.0;\n')
        f.write('    format      ascii;\n')
        f.write('    class       vectorField;\n')
        f.write('    location    "constant/polyMesh";\n')
        f.write('    object      points;\n')
        f.write('}\n')
        f.write('// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n')
        f.write('\n')
        f.write('\n')

        nodes = self.getCommonGrid()
        nodes = nodes.reshape((len(nodes)/3, 3))
        nPoints = len(nodes)
        f.write('%d\n'% nPoints)
        f.write('(\n')
        for i in range(nPoints):
            f.write('(%20.12f %20.12f %20.12f)\n'%(nodes[i, 0], nodes[i, 1], nodes[i, 2]))
        f.write(')\n\n\n')
        f.write('// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n')

    def printCurrentOptions(self):
        """Prints a nicely formatted dictionary of all the current options to
        the stdout on the root processor
        """
        if self.comm.rank == 0:
            print('+---------------------------------------+')
            print('|     All pyWarpUstruct Options:        |')
            print('+---------------------------------------+')
            pprint(self.solverOptions)

# =========================================================================
#                     Internal Private Functions
# =========================================================================
    def _setInternalSurface(self):
        """
        This function is used by default if setExternalSurface() is
        not set BEFORE an operation is requested that requires this
        information. Currently not written.
        """
        pass

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
        self.warp.gridinput.errtol = o['errTol']
        self.warp.kd_tree.bucket_size = o['bucketSize']
        if o['evalMode'].lower() == 'fast':
            self.warp.gridinput.evalmode = self.warp.gridinput.eval_fast
        elif o['evalMode'].lower() == 'exact':
            self.warp.gridinput.evalmode = self.warp.gridinput.eval_exact

    def __del__(self):
        """Release all the mesh warping memory"""
        self.warp.releasememory()

    def _getOFFileNames(self, caseDir):
        """
        Generate the standard set of filename for an openFoam grid.

        Parameters
        ----------
        caseDir : str
            The folder where the files are stored. 
        """
        if self.comm.size > 1:
            self.OFData['refPointsFile'] = os.path.join(caseDir, 'processor%d/constant/polyMesh/points_orig'%self.comm.rank)
            self.OFData['pointsFile'] = os.path.join(caseDir, 'processor%d/constant/polyMesh/points'%self.comm.rank)
            self.OFData['boundaryFile'] = os.path.join(caseDir, 'processor%d/constant/polyMesh/boundary'%self.comm.rank)
            self.OFData['faceFile'] = os.path.join(caseDir, 'processor%d/constant/polyMesh/faces'%self.comm.rank)
            self.OFData['ownerFile'] = os.path.join(caseDir, 'processor%d/constant/polyMesh/owner'%self.comm.rank)
            self.OFData['neighbourFile'] = os.path.join(caseDir, 'processor%d/constant/polyMesh/neighbour'%self.comm.rank)
        else:
            self.OFData['refPointsFile'] = os.path.join(caseDir, 'constant/polyMesh/points_orig')
            self.OFData['pointsFile'] = os.path.join(caseDir, 'constant/polyMesh/points')
            self.OFData['boundaryFile'] = os.path.join(caseDir, 'constant/polyMesh/boundary')
            self.OFData['faceFile'] = os.path.join(caseDir, 'constant/polyMesh/faces')
            self.OFData['ownerFile'] = os.path.join(caseDir, 'constant/polyMesh/owner')
            self.OFData['neighbourFile'] = os.path.join(caseDir, 'constant/polyMesh/neighbour')

    def _parseFoamHeader(self, lines, i):
        """Generic function to read the openfoam file header up to the point
        where it determines the number of subsequent 'stuff' to read
        
        Parameters
        ----------
        lines : list of strings
            Lines of file obtained from f.readlines()
        i : int
            The starting index to look
        
        Returns
        -------
        foamDict : dictionary
             Dictionary of the data contained in the header
        i : int 
             The next line to start at
        N : int
             The number of entries to read
        """
        keyword = re.compile(r'\s*[a-zA-Z]{1,100}\s*\n')

        blockOpen = False
        foamDict = {}
        imax = len(lines)
        while i < imax:
            res = keyword.match(lines[i])
            if res:
                i, foamData = self._readBlock(lines, i)
                break
            else:
                i += 1
        if i == imax:
            raise Error("Error parsing openFoam file header.")

        # Now we need to match the number followed by an open bracket:
        numberHeader = re.compile(r'\s*(\d{1,100})\s*\n')
        while i < imax:
            res = numberHeader.match(lines[i])
            if res:
                N = int(lines[i][res.start(0):res.end(0)])
                # We return +2 lines to skip the next '('
                return foamDict, i+2, N
            i += 1
        raise Error("Could not find the starting data in openFoam file")

    def _readBlock(self, lines, i):
        """
        Generic code to read an openFoam type block structure
        """

        openBracket = re.compile(r'\s*\{\s*\n')
        closeBracket = re.compile(r'\s*\}\s*\n')
        data = re.compile(r'\s*([a-zA-Z]*)\s*(.*);\s*\n')

        blockOpen = False
        imax = len(lines)
        blk = {}
        while i < imax:
            if not blockOpen:
                res = openBracket.match(lines[i])
                if res:
                    blockOpen = True
                i += 1
            else:
                res = data.match(lines[i])
                if res:
                    blk[res.group(1)] = res.group(2)
                if closeBracket.match(lines[i]):
                    i += 1
                    break
                i += 1
        return i, blk

    def _getCGNSFamilyList(self):
        """Obtain the family list from fortran. Only the root proc has this,
        so we need to bcast. 
        """
        familyList = None
        if self.comm.rank == 0:
            fullList = self.warp.griddata.familylist.transpose().flatten()
            nFamilies = self.warp.griddata.nwallfamilies
            familyList = ["" for i in range(nFamilies)]
            for i in range(nFamilies):
                tempName = fullList[i*32:(i+1)*32]
                newName = ''.join([tempName[j] for j in range(32)]).lower()
                familyList[i] = newName.strip()
        familyList = self.comm.bcast(familyList)

        # Converto back to a dict and assume all are walls:
        self.familyList = {}
        for fam in familyList:
            self.familyList[fam] = IWALL

    def _getIndices(self, groupName):
        """
        Try to see if groupName is already set, if not raise an
        exception

        """
        try: 
            indices = self.familyGroup[groupName]['indices']
        except:
            raise Error('%s has not been set using addFamilyGroup'%groupName)

        return indices

    def _checkOptions(self, solverOptions):
        """
        Check the solver options against the default ones
        and add opt
        ion iff it is NOT in solverOptions

        """
        for key in self.solverOptionsDefault.keys():
            if not key in solverOptions.keys():
                solverOptions[key] = self.solverOptionsDefault[key]
            else:
                self.solverOptionsDefault[key] = solverOptions[key]	

        return solverOptions
