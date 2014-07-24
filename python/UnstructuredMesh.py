#!/usr/bin/python
from __future__ import print_function
"""
Unstructure Mesh

The UnstructuredMesh module is used for interacting with an unstructured
mesh - typically used in a 3D CFD program.

It contains the following classes:

USMesh: General class for working with unstructured meshes

Copyright (c) 2014 by C.A.(Sandy) Mader
All rights reserved. Not to be used for commercial purposes.
Revision: 1.0   $Date: 09/07/2014$

Developers:
-----------
- C.A.(Sandy) Mader (CAM)

History
-------
	v. 1.0 - Initial Class Creation (CAM, 2014)

"""
# =============================================================================
# Standard Python modules
# =============================================================================
import sys, os, time
import subprocess,re

# =============================================================================
# External Python modules
# =============================================================================
import numpy as np
try:
    # import the necessary modules from PyFoam
    from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
    from PyFoam.RunDictionary.ListFile import ListFile
    from PyFoam.RunDictionary.MeshInformation import MeshInformation
except:
    print('PyFoam is unavailable, OpenFOAM case will not run...')
# end


# =============================================================================
# Extension modules
# =============================================================================
from pygeo import geo_utils
from mpi4py import MPI
from MExt import MExt

class Error(Exception):
    """
    Format the error message in a box to make it clear this
    was a expliclty raised exception.
    """
    def __init__(self, message):
        msg = '\n+'+'-'*78+'+'+'\n' + '| pyWarpUstruct Error: '
        i = 23
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

class USMesh(object):
    """
    This is the main Unstructured Mesh. This mesh object is designed to
    interact with an unstructured CFD solver though a variety of
    interface functions.
 
    """
    def __init__(self, fileName, comm=None, meshOptions=None,*args,**kwargs):
        """
        Create the USMesh object. 

        Input Arguments:
            fileName, str: either a CGNS file name or and OpenFOAM case
                           directory from which to create the mesh
        Optional Arguments:
            comm, MPI_INTRA_COMM: MPI communication (as obtained 
                from mpi4py) on which to create the USMesh object.
                If not provided, MPI_COMM_WORLD is used. 
             meshOptions, dictionary: A dictionary containing the 
                 the options for the mesh movement strategies. 

         Returns:
             None

             """

        if comm is None:
            self.comm = MPI.COMM_WORLD
        else:
            self.comm = comm
            
        # Check if warp has already been set. If this is true, this has been
        # inherited to the complex version
        try: 
            self.warp
        except:
            if 'debug' in kwargs:
                debug = kwargs['debug']
            else:
                debug = False
            curDir = os.path.dirname(os.path.realpath(__file__))
            self.warp = MExt('warpustruct',path=[curDir], debug=debug)._module
        # end try

        if meshOptions is not None:
            self.solverOptions = meshOptions
        else:
            self.solverOptions = {}
        # end if

        # Set realtype of 'd'. 'D' is used in Complex and set in
        # UnstructuredMesh_C.py
        self.dtype = 'd'

        # Defalut options for mesh warping
        self.solverOptionsDefault = {
            # Warp Type
            'warpType':'inversedistance',#alt: rbf
            'aExp': 3,
            'bExp': 5,
            'fileType':'cgns', # other option openfoam
            'binaryPointsFile':False,
            'binaryFacesFile':False,
            'binaryCellsFile':False,
            'binaryOwnersFile':False,
            }
        
        # # Set remaining default values
        # self.FETopo = None
        self._checkOptions(self.solverOptions)
        #self._setMeshOptions()
        # self.warpMeshDOF = self.warp.block.nnodeslocal*3
        # self.solverMeshDOF = 0
        self.warpInitialized = False
        self.IDWarpInitialized = False


        # Read the grid from the CGNS file
        print(fileName)
        fileName, ext = os.path.splitext(fileName)
        print('filename after',fileName,ext)
       
        if self.solverOptions['fileType'] == 'cgns':
            self.warp.readunstructuredcgnsfile(fileName+'.cgns', self.comm.py2f())

        elif self.solverOptions['fileType'] == 'openfoam':
            # we are dealing with an OpenFOAM case
            self._readOpenFOAMGrid(fileName)
        else:
            print('Invalid input file type. CGNS or OpenFOAM required...Exiting')
            sys.exit(0)
        # end if
       
        # Figure out the unique surface numbering
        self.preprocessSurfaces()

        # Extract family info and add default group "all"
        self._getFamilyList()
        self.familyGroup = {}
        print('adding family group')
        self.addFamilyGroup('all')


        return

# =========================================================================
#                  Local mesh reading functionality
# =========================================================================

    def _readOpenFOAMGrid(self,caseDir):
        '''
        Read in the mesh points and connectivity from the polyMesh directory
        '''
        
        # Copy the reference points file to points to ensure consistant starting
        # point
        self.refPointsFile = os.path.join(caseDir,'constant/polyMesh/points_orig')
        self.pointsFile = os.path.join(caseDir,'constant/polyMesh/points')

        status = subprocess.call("cp " + "%s %s"%(self.refPointsFile,self.pointsFile), shell=True)
        if status != 0:
            raise Error('pyWarpUstruct: status %d: Unable to copy points_orig to points.'%status)

        # Create an instance of mesh info
        self.mInfo = MeshInformation(caseDir)

        # Read in the volume points
        self.nodes = self._readVolumeMeshPoints()

        # Read the boundary info
        self._readBoundaryInfo()
        
        # Read the face info for the mesh
        self._readFaceInfo()

        # Read the cell info for the mesh
        self._readCellInfo()

        return

    def _readVolumeMeshPoints(self):
        '''
        return an numpy array of the mesh points
        '''

        # use PyFoam to determine the number of points in the file
        self.nPoints = self.mInfo.nrOfPoints()
        #print 'nPoints',self.nPoints
       
        # Open the points file for reading
        pointHandle = open(self.pointsFile,'r')

        self._parseHeaderNumerical(pointHandle,self.nPoints)
      
        # now read the points into x
        print (self.solverOptions)
        print ('binaryPointsFile?',self.solverOptions['binaryPointsFile'])
        if self.solverOptions['binaryPointsFile']:
            # read the points in binary using fromfile
            x = np.fromfile(pointHandle, dtype='float', count=self.nPoints*3, sep="").reshape((self.nPoints,3))

        else:
            # read the file in as ascii
            x = np.zeros([self.nPoints,3])
            counter = 0
            for line in pointHandle.readlines():
                #print 'before',line
                line = re.sub('[)(]', '', line)
                #print 'after',line
                vals = line.split()
                #print len(vals),vals
                if len(vals)==3 and counter<self.nPoints:
                    #print 'points found',counter
                    x[counter,0] = float(vals[0])
                    x[counter,1] = float(vals[1])
                    x[counter,2] = float(vals[2])
                    counter += 1
                # end
            # end
        # end
                    
        pointHandle.close()

        return x

    def _readBoundaryInfo(self):
        '''
        read the boundary file information for this case and store in a dict.
        '''
        dirName = os.getcwd()
        
        boundaryFile = os.path.join(dirName,'constant/polyMesh/boundary')
        meshDir = os.path.join(dirName,'constant/polyMesh/')

        boundary=ListFile(meshDir,"boundary")
        
        nBoundaries = boundary.getSize()

        boundaryHandle = open(boundaryFile, 'r')

        self._parseHeaderNumerical(boundaryHandle,nBoundaries)

        print ('boundary',dir(boundary),boundary.getSize())
        
        boundaries = {}
        boundaryKeys = []
        boundaryCounter = 0
        startDict = True
        for line in boundaryHandle.readlines():
            #print line
            vals = line.split()
            if len(vals)==0:
                continue
            else:
                if startDict:
                    boundaries[vals[0]] = {}
                    boundaryKeys.append(vals[0])
                   
                    startDict = False
                elif boundaryCounter<nBoundaries-1:
                    if '{' in line:
                        continue
                    elif '}' in line:
                        startDict=True
                        boundaryCounter+=1
                        print(boundaryCounter,nBoundaries)
                    else:
                        boundaries[boundaryKeys[boundaryCounter]][vals[0]] = vals[1]
                    # end
                # end
                #print vals

            # end
        # end
        print(boundaries.keys())

        boundaryHandle.close()
        for key in boundaries.keys():
            print(boundaries[key])
                
        sys.exit(0)
        return

    def _readFaceInfo(self):
        '''
        Read the face info for this case.
        '''
        dirName = os.getcwd()
        
        faceFile = os.path.join(dirName,'constant/polyMesh/faces')
        #faceFile2 = os.path.join(dirName,'constant/polyMesh/faces_test')
        meshDir = os.path.join(dirName,'constant/polyMesh/')

        faces=ListFile(meshDir,"faces")

        self.nFaces = faces.getSize()

        faceHandle = open(faceFile, 'r')
        #faceHandle2 = open(faceFile2, 'rb')

        # parse through the header to get to the start of the binary data
        self._parseHeaderNumerical(faceHandle,self.nFaces)

        if self.solverOptions['binaryFacesFile']:
            for  i in range(30):
                print(ord(faceHandle2.read(1)))
            # end

            for i in range(2):
                value = np.fromfile(faceHandle2, dtype='short', count=1, sep="")
                print('val',value)
                # values = np.fromfile(faceHandle2, dtype='<i4', count=value, sep="")
                # print values
            sys.exit(0)
            for i in range(80):#5):#nFaces):
                print(ord(faceHandle2.read(1)))

            # for i in range(50):
            #     #value = faceHandle2.read(8)
            #     value = np.fromfile(faceHandle, dtype='int', count=1, sep="")
            #     print 'val',value
            
            # for i in range(2):#5):#nFaces):
            #     print faceHandle.read(1)
            # for i in range(6):#5):#nFaces):
            #     nFacePoints = np.fromfile(faceHandle, dtype='int', count=nFaces, sep="")
            #     print nFacePoints
            #     # print nFacePoints
                
            #     # pointIndices = np.fromfile(faceHandle, dtype=np.int32, count=nFacePoints, sep="")
            #     # print pointIndices


            # end
            sys.exit(0)
        else:
            # read the file in as ascii
            faces = {}
            counter = 0
            for line in faceHandle.readlines():
                line = re.sub('[)(]', ' ', line)
                vals = line.split()
                if len(vals)==0:
                    continue
                elif "//" in line:
                    continue
                elif counter<self.nFaces:
                    nPoints = int(vals[0])
                    faces[counter]= np.zeros(nPoints,int)
                    for i in range(nPoints):
                        faces[counter][i] = int(vals[i+1])
                    # end
                    counter += 1
                # end
            # end

        self.faces = faces
        #print 'faces',dir(faces),faces.getSize()
        
        # tmp = faces.readFile()
        # print tmp
        # self.faces=faces.getSize()

        faceHandle.close()
        return

    def _readCellInfo(self):
        '''
        read the boundary file information for this case and store in a dict.
        '''
        dirName = os.getcwd()
        
        cellFile = os.path.join(dirName,'constant/polyMesh/cells')
        ownerFile = os.path.join(dirName,'constant/polyMesh/owner')
        meshDir = os.path.join(dirName,'constant/polyMesh/')

        cells=ListFile(meshDir,"cells")
        owners = ListFile(meshDir,"owner")
        
        self.nOwners = owners.getSize()
        self.nCells = cells.getSize()

        cellHandle = open(cellFile, 'r')
        ownerHandle = open(ownerFile, 'r')

        # parse through the header to get to the start of the binary data
        self._parseHeaderNumerical(cellHandle,self.nCells)
        self._parseHeaderNumerical(ownerHandle,self.nOwners)

        # read the cells
        if self.solverOptions['binaryCellsFile']:
            pass
        else:
            # read the file in as ascii
            cells = {}
            counter = 0
            for line in cellHandle.readlines():
                line = re.sub('[)(]', ' ', line)
                vals = line.split()
                if len(vals)==0:
                    continue
                elif "//" in line:
                    continue
                elif counter<self.nCells:
                    nFacesPerCell = int(vals[0])
                    cells[counter]= np.zeros(nFacesPerCell,int)
                    for i in range(nFacesPerCell):
                        cells[counter][i] = int(vals[i+1])
                    # end
                    counter += 1
                # end
            # end

        # read the owners
        if self.solverOptions['binaryOwnersFile']:
            pass
        else:
            # read the file in as ascii
            owners = {}
            self.owners = np.zeros(self.nOwners,int)
            counter = 0
            for line in ownerHandle.readlines():
                line = re.sub('[)(]', ' ', line)
                vals = line.split()
                if len(vals)==0:
                    continue
                elif "//" in line:
                    continue
                elif counter<self.nOwners:
                    self.owners[counter]=int(vals[0])
                    counter += 1
                # end
            # end

        

    def _parseHeaderNumerical(self,handle,stopValue):
        '''
        a generic function to parse through a file and return when a certain
        number is detected in the file
        '''

        # Setup a character based check to find the start of the list entries
        # Turn the number of entries into a string
        nEntriesChar = '%d'%stopValue
        # determine the length of the string
        lenNEntries = len(nEntriesChar)
        # create a boolean array with an entry for each character of the string
        # if we match all of the characters consecutively, we will set a logical
        # which will combine with the '(' character to indicate the start of the list
        entriesCheck = np.zeros(lenNEntries,bool)

        # set and index for the character boolean and create the logical to indicate
        # when we have found the consecutive string of characters that represent the
        # number of points
        eIdx = 0
        nEntriesFound = False
        # loop over the characters in the file one byte at a time
        while True:
            c = handle.read(1)
            if not nEntriesFound:
                # check to see if the current character matches the required 
                # character in the sequence
                if c == nEntriesChar[eIdx]:
                    # if they match, set the current entry in the boolean array 
                    # to true and increment the index to the next entry in the 
                    # array
                    print(c)
                    entriesCheck[eIdx] = True
                    eIdx += 1

                else:
                    # Otherwise we haven't found the string so we need to start over
                    # reset everything
                    entriesCheck[:] = False
                    eIdx = 0

                # end

                # check if all of the entries in the boolean array are True
                # if so we have found the nEntries string and we can start searching
                # for the start of the list of entries
                if all(entriesCheck):
                    nEntriesFound = True
                # end
                        
            else:
                # Now that we have found the nEntries Character string, we can start
                # looking for the '(' that represents the start of the entries list
                # once we have found that character, exit the loop.
                if c=='(':
                    print("final Character",c)
                    print("End of header")
                    break
                # end
            # end
            #print "Read a character:", c
        # end

        return
        
# ============================
# Local grid geometry calcs
# ============================
    def preprocessSurfaces(self):
        '''
        Run the routines to separate out the surfaces and determine the 
        node orderings
        '''

        #self.warp.getuniquesurfacenodelist()
        self.warp.getfulluniquesurfacenodelist()
        self.updateGridMetrics(True)

        return

    def updateGridMetrics(self,initialPoint=False):
        '''
        run the fortran level routines to compute area for 
        the surface elements/faces
        '''

        self.warp.getsurfaceelementcenterandarea()
        self.warp.computenodalproperties(initialPoint)


        return

# =========================================================================
#                  Local Multi-disciplinary Surface Functionality
# =========================================================================

    def getSurfaceCoordinates(self, groupName):
        """ 
        Returns all surface coordinates on this processor in group
        'groupName'

        Input Arguments:
           groupName, str: The group from which to obtain the coordinates.
           This name must have been obtained from addFamilyGroup() or
           be the default 'all' which contains all surface coordiantes
        Output Arguements:
            coords, numpy array, size(N,3): coordinates of the requested
            group. This may be empty arry, size (0,3)
            
        """
        self._setInternalSurface()
        indices = self._getIndices(groupName)
        coords = np.zeros((len(indices),3), self.dtype)
        print('about to get coords',groupName,len(indices))
        self.warp.getsurfacecoordinates(indices, np.ravel(coords))
        print('points',coords.shape)

        return coords
   
#     def getSurfaceConnectivity(self, groupName):
#         """
#         Returns the connectivities of the surface coordinates that
#          were obtained with getSurfaceCoordinates()
#         """

#         return self.familyGroup[groupName]['connectivity']

    def setSurfaceCoordinates(self, groupName, coordinates):
        """ Sets all surface coordinates on this processor in group
        "groupName"
        
        Input Arguments:
           groupName, str: The group to set the coordinate in. 
               This name must have been obtained from addFamilyGroup() or
               be the default 'all' which contains all surface coordiantes
           coordinates, numpy array, size(N,3): The coordinate to set. This MUST
               be exactly the same size as the array obtained from 
               getSurfaceCoordinates()
         Output Arguements:
            Noneupd

            """

        indices = self._getIndices(groupName)
        self.warp.setsurfacecoordinates(indices, np.ravel(coordinates))
        
        return 
 
    def addFamilyGroup(self, groupName, families=None):
        """ 
        Create a grouping of CGNS families called "groupName"
        
        Input Arguments:
            groupName, str: The name to call this collection of families
        Optional Arguments:
            families, list: A list containing the names of the CGNS familes
                the user wants included in "groupName". The items in the list
                must be valid family names. The user can call printFamilyList() 
                to determine what families are present in the CGNS grid.
                Default: All families
        Output Arguements:
            None

            """
      
        # Use whole list by default
        if families == None: 
            families = self.familyList
        # end if

        groupFamList = []
        for fam in families:
            if fam.lower() in self.familyList:
                groupFamList.append(fam)
            else:
                print('* WARNING: %s was not in family list. \
I will ignore this family'%(fam),comm=self.comm)
            # end if
        # end for
        
        self.familyGroup[groupName] = {'families':groupFamList}
#        print ('family Group',groupName,self.familyGroup[groupName])
        if self.warpInitialized: # We can add the indices
            nodeIndices = []
            nodeCount = 0

            for i in xrange(len(self.patchNames)):
                sizeOfPatch = self.patchSizes[i]#[0]*self.patchSizes[i][1]*3
                
                if self.patchNames[i] in groupFamList:
                    # Get the nodes indices:
                    nodeIndices.extend(np.array(self.patchIndices[i]))
#                    print('nodeLength',len(self.patchIndices[i]),len(nodeIndices))
                # end if
                nodeCount += sizeOfPatch
            # end for
#            print('groupName',groupName,nodeCount,len(nodeIndices))
            self.familyGroup[groupName]['indices'] = nodeIndices
            
        #     cellConn = np.zeros((0,4),'intc')
        #     if self.connectivity is None:
        #         self.familyGroup[groupName]['connectivity'] = cellConn
        #         return 
         
        #     nNodesSkipped = 0
        #     cCnt = 0
          
        #     for i in xrange(len(self.patchNames)):
        #         cellSize = (self.patchSizes[i][0]-1)*(self.patchSizes[i][1]-1)            
        #         if self.patchNames[i] in groupFamList: 
        #             # Add the connecitivity minus the number of nodes
        #             # we've skipped so far. 
        #             tmp = self.connectivity[cCnt:cCnt+cellSize,:].copy() - nNodesSkipped
        #             cellConn = np.append(cellConn, tmp)
        #         else:
        #             # If we didn't add this patch increment the number
        #             # of nodes we've skipped
        #             nNodesSkipped += self.patchSizes[i][0]*self.patchSizes[i][1]
        #         # end if
        #         cCnt += + cellSize
        #     # end for

          
        #     self.familyGroup[groupName]['connectivity'] = cellConn

        return
   
    def printFamilyList(self):
        """ Prints the families in the CGNS file"""
        
        if self.comm.rank == 0:
            print('Family list is:', self.familyList)

        return

# # =========================================================================
# #                         Interface Functionality
# # =========================================================================

#     def setExternalMeshIndices(self, ind):
#         """ 
#         Set the indicies defining the transformation of an external
#         solver grid to the original CGNS grid. This is required to use
#         MBMesh functions that involve the word "Solver" and warpDeriv
        
#         Input Arguments:
#             ind, numpy integer array: The list of indicies this processor
#                 needs from the CGNS mesh

#                 """
#         self.warp.setexternalmeshindices(ind)
#         self.solverMeshDOF = len(ind)

#         return 

    def setExternalSurface(self, patchNames, patchSizes, patchIndices, conn, pts):
        """
        This is the Master function that defines the surface
        coordinates the external solver wants to use. The actual
        coordinates themselves MUST match exactly with the warp's
        surface (i.e. you must use the same mesh), but they need not
        match by processor or and the solver may have split pataches. 
        """

        self.patchNames = patchNames
        self.patchSizes = patchSizes
        self.patchIndices = patchIndices
        self.connectivity = conn
        
        # Call the fortran initialze warping command withe the
        # coordiantes we now have'

        # self.warp.initializewarping(np.ravel(pts.real.astype('d')))
        # self._initializeSolidWarping()
        self.warpInitialized = True

        # We can now back out the indices that should go along with
        # the groupNames that may have already been added. 
        print('running add families')
        for key in self.familyGroup.keys():
            self.addFamilyGroup(key, self.familyGroup[key]['families'])
        # end for
        
        return 

    def _setInternalSurface(self):
        """
        This function is used by default if setExternalSurface() is not
        set BEFORE an operation is requested that requires this information. 
        """

        if self.warpInitialized is False:
            print(' -> Info: Using Internal pyWarp Surfaces')
            # get the number of wall patches on this proc
            npatch = self.warp.getnpatches()
            patchNames = []
            patchIndices = []
            patchSectionIndices = []
            patchSizes = []
            ptSize = 0
            #print ('Npatches',npatch)
            # loop over the patches to get their names and sizes
            for i in xrange(npatch):               
                tmp = self.warp.getpatchname(i)
                patchIdx = self.warp.getpatchindex(i)
                #print ('patchname',i,tmp,patchIdx)
                patchSectionIndices.append(patchIdx)
                patchNames.append(
                    ''.join([tmp[j] for j in range(32)]).lower().strip())
                patchSizes.append(self.warp.getpatchsize(patchSectionIndices[i]))
                #print('patchSizes',patchSizes[i])
                patchIndices.append(self.warp.getpatchindexlist(patchSectionIndices[i],patchSizes[i]))
                #print('index length',len(patchIndices[i]))
                ptSize += patchSizes[-1]#[0]
            # end for

            # get the patch points
            conn = None
            pts = np.zeros((ptSize,3), 'd') #Explicitly real...even in complex
            self.warp.getpatches(np.ravel(pts))
            #print(pts.shape,patchNames,patchSizes[0])

            # Run the "external" command
            self.setExternalSurface(patchNames, patchSizes, patchIndices, conn, pts)

            self.warpInitialized = True
        return

#     def getSolverGrid(self):
#         """
#         Return the current grid in the order specified by
#         setExternalMeshIndices()

#         Input Arguments: 
#             None
#         Output Arguments:
#             solverGrid, numpy array, real: The resulting grid. 
#                 The output is returned in flatted 1D coordinate
#                 format. The len of the array is 3*len(indices) as
#                 set by setExternalMesIndices()
#                 """

#         solverGrid = np.zeros(self.solverMeshDOF, self.dtype)
#         warpGrid   = np.zeros(self.warpMeshDOF, self.dtype)
#         self.warp.packblocks(warpGrid)
#         self.warp.warp_to_solver_grid(warpGrid, solverGrid)
        
#         return solverGrid

#     def getWarpGrid(self):
#         """
#         Return the current grid. This funtion is typically unused. See
#         getSolverGrid for the more useful interface functionality. 

#         Input Arguments: 
#             None
#         Output Arguments:
#             warpGrid, numpy array, real: The resulting grid. 
#                 The output is returned in flatted 1D coordinate
#                 format. 
#                 """

#         warpGrid = np.zeros(self.warpMeshDOF, self.dtype)
#         self.warp.packblocks(warpGrid)
        
#         return warpGrid

#     def getdXs(self, groupName):
#         """
#         Return the current values in dXs
        
#         Input Arguments:
#             groupName, str: The group defining which components will be
#                returned from dXs
#         Output Arguments:
#             dXs, numpy array real: The specific components of dXs. 
#                 size(N,3) where N is the number of points associated
#                 with the families defined by "groupName"
                
#                 """

#         indices = self._getIndices(groupName)
#         dXs = np.zeros((len(indices)/3,3),self.dtype)
#         self.warp.getdxs(indices, np.ravel(dXs))

#         return dXs

#     def setdXs(self, dXs, groupName):
#         """
#         Return the current values in dXs
        
#         Input Arguments:
#             dXs, real numpy array: The components correponding to
#                 groupName to set in dXs. 
#             groupName, str: The group defining where components will be
#                 set in dXs
#         Output Arguments:
#             None
            
#             """
        
#         indices = self._getIndices(groupName)
#         self.warpsetdxs(indices, np.ravel(dXs))

#         return

#     def sectionVectorByFamily(self, groupName, inVec):
#         """
#         Partition a dXs like vector (inVec) that is "full size", ie
#         size of the "all" group and return just the part corresponding
#         to "groupName"
#         Input Arguments:
#             groupName, str: The group defining where components will be
#                 returned
#             inVec, real numpy array: The components to section whose length
#                 corresponds to the len of the "all" group
#                 groupName to set in dXs. 
#         Output Arguments:
#             dXs: The sectioned form of inVec
            
#             """
#         indices = self._getIndices('all')
#         self.warp.setdxs(indices, np.ravel(inVec))
#         indices = self._getIndices(groupName)
#         outVec = np.zeros((len(indices)/3,3), self.dtype)
#         self.warp.getdxs(indices, np.ravel(outVec))

#         return outVec

#     def expandVectorByFamily(self, groupName, inVec):
#         """
#         Expand a dXs like vector (inVec) that correponds to
#         "groupName" and return a vector that is of full size, ie the
#         "all" group. It does the reverse operation of
#         sectionVectorByFamily(). All other entries of the output are
#         zeroed
#         Input Arguments:
#             groupName, str: The group defining where components will be
#                 added
#             inVec, real numpy array: The components to section whose length
#                 corresponds to the len of "groupName"
#         Output Arguments:
#             dXs: The expanded form of dXs corresponding to the "all"
#                 group
                
#                 """
#         indices = self._getIndices(groupName)
#         self.warp.setdxs(indices, np.ravel(inVec))
#         indices = self._getIndices('all')
#         outVec = np.zeros((len(indices)/3,3), self.dtype)
#         self.warp.getdxs(indices, np.ravel(outVec))
        
#         return outVec

# # ==========================================================================
# #                        Output Functionality
# # ==========================================================================
    def writeGridTecplot(self,fileName):
        '''
        write the current grid coordinates in a tecplot FE
        File.
        '''
        f = open(fileName+'.dat','w')
        
        nPoints = self.nPoints
        nFaces = self.nFaces
        nCells = self.nCells

        # write the node numbers for each face
        faceNodeSum = 0
        for i in range(nFaces):   
            nPointsFace = len(self.faces[i])
            faceNodeSum+=nPointsFace
        # end

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
        for i in range(3):
            for j in range(nPoints):
                f.write('%f\n'%self.nodes[j,i])
            # end
        # end
        # Write the number of nodes in each face to the
        # file
        counter = 0
        for i in range(nFaces):
            nPointsFace = len(self.faces[i])
            f.write('%d '%nPointsFace)
            counter += 1
            if counter>300:
                f.write('\n')
                counter=0
        # end
        f.write('\n')
        # write the node numbers for each face
        for i in range(nFaces):   
            nPointsFace = len(self.faces[i])
            for j in range(nPointsFace):
                f.write('%d '%(self.faces[i][j]+1))
            # end
            f.write('\n')
        # end
        # write left elements
        for i in range(nFaces):
            f.write('%d\n'%0)#self.owners[i])
        # write right elements
        for i in range(nFaces):
            f.write('%d\n'%1)#self.owners[nFaces-i-1])
            

        f.close()
        return
#     def writeVolumeGrid(self, fileName):
#         """
#         Write the current state the mesh to a volume CGNS file
        
#         Input Arguments: 
#             fileName: cgns file name for output file

#         Output Arguments:
#             None

#             """
        
#         self.warp.writevolumegrid(fileName)

#         return

#     def writeSuperNodeVolumeGrid(self, fileName):
#         """
#         Write the current state the super node mesh to a volume CGNS
#         file
        
#         Input Arguments: 
#             fileName: cgns file name for output file

#         Output Arguments:
#             None

#             """ 

#         self._setInternalSurface()
#         if self.solverOptions['warpType'] != 'algebraic':
#             self.warp.writevolumegrid_super(fileName)
#         # end if

#         return

#     def writeXswVolumeGrid(self, fileName):
#         """
#         Write the solid warp fill approximation of the volume grid. This is
#         normally only used for illustrative purposes. 

        
#         Input Arguments: 
#             fileName: cgns file name for output file

#         Output Arguments:
#             None

#             """ 

#         self._setInternalSurface()
#         if self.solverOptions['warpType'] != 'algebraic':
#             self.warp.swap_x_and_xsw()
#             self.warp.writevolumegrid(fileName)
#             self.warp.swap_x_and_xsw()
#         # end if
        
#         return

#     def writeSurfaceGrid(self, fileName):
#         """
#         Write the current state the surface mesh to a volume CGNS file
        
#         Input Arguments: 
#             fileName: cgns file name for output file

#         Output Arguments:
#             None

#             """

#         self.warp.writesurfacegrid(fileName)

#         return

#     def writeSuperNodeSurfaceGrid(self, fileName):
#         """
#         Write the current state the Super Node surface mesh to a
#         surface CGNS file
        
#         Input Arguments: 
#             fileName: cgns file name for output file

#         Output Arguments:
#             None

#             """ 

#         self._setInternalSurface()
#         if self.solverOptions['warpType'] != 'algebraic':
#             self.warp.writesurfacegrid_super(fileName)
#         # end if

#         return

#     def writeXswSurfaceGrid(self, fileName):
#         """
#         Write the solid warp fill approximation of the surface
#         grid. This is normally only used for illustrative purposes.

        
#         Input Arguments: 
#             fileName: cgns file name for output file

#         Output Arguments:
#             None

#             """ 

#         self._setInternalSurface()
#         if self.solverOptions['warpType'] != 'algebraic':
#             self.warp.swap_x_and_xsw()
#             self.warp.writesurfacegrid(fileName)
#             self.warp.swap_x_and_xsw()
#         # end if

#         return
        
#     def writeTopology(self, fileName):
#         """
#         Write the computed grid topology to a file

#         Input Arguments: 
#             fileName, str: The file name for the .con file

#         Usage Notes: This function write out the grid connectivity as
#         well an ASCII Tecplot file defining the edges. The idea is to
#         run your script once, use this function to create the topo
#         file along with the design group labels. The design group
#         labels can be loaded into tecplot along with the volume grid
#         and the user can determine the desired number of "Super Nodes"
#         along each design group. The desired numbers can then be
#         entered in the topology file.  On subsequent runs, the user
#         uses the options:

#         'solidWarpType':'topo', 
#         'topo':'topo_file.con'

#         to use the supplied supplied topology file.
        
#         """
#         if self.FE_topo is not None:
#             if self.comm.rank == 0:
#                 self.FETopo.writeConnectivity(fileName + '.con')
            
#             # Get the set of corners and mid points from the file:
#             nBlocks, coords, BCs, blockDims = self._getCGNSCoords()        

#             # Next write out edge file:
#             if self.comm.rank == 0:
#                 labelFilename = fileName + '.dat'
#                 f = open(labelFilename, 'w')
#                 for iVol in xrange(nBlocks):
#                     # We ONLY want to write the edge if its on a
#                     # surface so as to not clutter the output.
#                     for iFace in xrange(6):
#                         if BCs[iVol, iFace] <> 0:
#                             edges = geo_utils.edgesFromFace(iFace)
#                             for iEdge in edges:
                        
#                                 ue = self.FETopo.edge_link[iVol, iEdge]
#                                 dg = self.FETopo.edges[ue].dg
#                                 pt = coords[iVol, 8+iEdge]
#                                 textString = 'TEXT CS=GRID3D X=%f,Y=%f,Z=%f,\
# T=\"DG%d\"\n'% (pt[0], pt[1], pt[2], dg)
#                                 f.write('%s'%(textString))
#                             # end for
#                         # end if
#                     # end for
#                 # end for
#                 f.close()
#             # end if
#         # end if

#         return 

# # =========================================================================
# #                      Mesh Warping Functionality
# # =========================================================================

    def _setMeshOptions(self):
        """ Private function to set the options currently in
        self.solverOptions to the corresponding place in Fortran"""
        
#         # Solid/Algebraic Warping
#         if self.solverOptions['warpType'] == 'algebraic':
#             self.warp.solidwarpmodule.warptype = \
#                 self.warp.solidwarpmodule.algebraicwarping
#         elif self.solverOptions['warpType'] == 'solid':
#             self.warp.solidwarpmodule.warptype = \
#                 self.warp.solidwarpmodule.solidwarping
#         else:
#             mpiPrint('warpType is not reconginzed. Applicable\
# values are \'solid\' and \'algebraic\'')
#             sys.exit(1)
#         # end if
            

        return

    def _initializeIDWarping(self):
        """
        Internal function to setup the necessary data for the inverse distance warping
        """
        
        if not self.IDWarpInitialized and \
                self.solverOptions['warpType'] == 'inverseDistance':

            print('\nInitializating InverseDistance Mesh Warping...',
                     comm=self.comm)

            # Determine the Boundary nodes

            # compute Ldef
            
            # Compute Alpha
            
            # Compute the area weighting of the boundary nodes
                         

            print('  -> Inverse Distance Mesh Warping Initialized.', comm=self.comm)
            self.IDWarpInitialized = True
        # end if

        return

    def warpMesh(self):
        """ 
        Run the applicable mesh warping strategy.

        This will update the volume coordinates to match surface
        coordinates set with setSurfaceCoordinates()
        
        """

        # Initialize internal surface if one isn't set
        self._setInternalSurface()

        # Set Options
        self._setMeshOptions()

        # Make sure solid warping is initialized if necessary
        self._initializeIDWarping()

        # Warp the mesh

        # Compute the updated surface element centers, areas and normals
        print('updating grid metrics')
        self.updateGridMetrics()

        # loop over the unique surface nodes to compute the rotations and
        # displacements
        
        # compute Si for eache surface node

        # compute Smean

        # Compute alpha for current displacement

        # For each volume node not on the boundary
        # loop over the surface and compute wi and si
        # Sum s(r) on the fly
        print('updateingVolume')
        self.warp.updatevolumecoordinates()
        
        return

#     def resetWarping(self):
#         """
#         This function is used to reset the non-linear warping
#         solution. If a bad solution is found and the mesh is invalid
#         this will reset the solution vector such that on a subsequent
#         movement a new starting point will be found using the stepped
#         approach with 'nSolutionSteps'
#         """
        
#         self.warp.resetwarping()

#         return

#     def warpDeriv(self, dXv, solverVec=True, surfOnly=False):
#         """
#         Compute the warping derivative (dXvdXs^T)*Vec

#         This is the main routine to compute the mesh warping
#         derivative. 

#         Input Arguments:

#             solverdXv, numpy array: Vector of size external
#                 solver_grid. This is typically obtained from the
#                 external solver's dRdx^T * psi calculation. 
        
#             solverVec, logical: Flag to indicate that the dXv vector
#                  is in the solver ordering and must be converted
#                  to the warp ordering first. This is the usual approach
#                  and thus defaults to True.

#             surfOnly, logical: Flag to indicate that a "fake" mesh
#                  warp is to be done. In this case, only dXv values that
#                  are ALREADY on the surface are taken. This is only used
#                  in a few select cases and thus defaults to False.

#         Output Arguments:

#             None. However, the resulting calculation is available from
#             the getdXs() function. 

#         """
#         if solverVec:
#             dXvWarp = np.zeros(self.warpMeshDOF, self.dtype)
#             self.warp.solver_to_warp_grid(dXv, dXvWarp)
#         else:
#             dXvWarp = dXv
#         # end if
#         if surfOnly:
#             self.warp.warpderivsurfonly(dXvWarp)
#         else:
#             self.warp.warpderiv(dXvWarp)
#         # end if

#         return 

#     def verifyWarpDeriv(self, dXv=None, solverVec=True, dofStart=0, 
#                         dofEnd=25):
#         """Run an internal verification of the solid warping
#         derivatives"""

        
#         if dXv is None:
#             np.random.seed(314) # 'Random' seed to ensure runs are same
#             dXvWarp = np.random.random(self.warpMeshDOF)
#         else:
#             if solverVec:
#                 dXvWarp = np.zeros(self.warpMeshDOF, self.dtype)
#                 self.warp.solver_to_warp_grid(dXv, dXvWarp)
#             else:
#                 dXvWarp = dXv
#         # end if

#         self.warp.verifywarpderiv(dXvWarp, dofStart, dofEnd)

#         return

#     def verifyWarpDerivFwd(self, callMode):
#         """Run an tapenade forward mode AD verification. Call mode
#         must be 1 for first call and 2 after"""

#         self.warp.verifywarpderivfwd(callMode)

#         return

#     def verifyWarpDerivBwd(self, callMode):
#         """Run an tapenade forward mode AD verification. Call mode
#         must be 1 for first call and 2 after"""

#         self.warp.verifywarpderivbwd(callMode)

#         return

# # =========================================================================
# #                         Utiliy Functions
# # =========================================================================

#     def getMeshQuality(self, bins=None):
#         """
#         Determine the mesh quality of the current volume mesh
        
#         Optional Arguments:
#             bins, numpy array: A list of numbers from -1 to -1 used 
#                 in binning the histogram. 

#         Output Arguments:
#             qualityMin, real: The minimum element quality
#             qualityMax, real: The maximum element quality
#             qualityAvg, real: The average element quality
#             hist, numpy array: The values for the resulting histogram
#             binEdges, numpy array: Values for the edges of the histogram

#             """
#         if bins is None:
#             bins = np.linspace(-1, 1, 21)
#         # end if

#         nvolLocal = self.warp.getnvolproc()
#         qualityLocal = self.warp.getquality(nvolLocal)
#         res = self.comm.gather(qualityLocal, root=0)

#         if self.comm.rank == 0:
#             quality = []
#             for i in range(self.comm.size):
#                 quality.extend(res[i])
#             # end for
                
#             qualityMin = np.min(quality)
#             qualityMax = np.max(quality)
#             qualityAvg = np.average(quality)

#             hist, binEdges = np.histogram(quality, bins=bins)
#         else:
#             qualityMin = None
#             qualityMax = None
#             qualityAvg = None
#             hist       = None
#             binEdges   = None
#         # end if
        
#         # Broadcast:
#         qualityMin = self.comm.bcast(qualityMin, root=0)
#         qualityMax = self.comm.bcast(qualityMax, root=0)
#         qualityAvg = self.comm.bcast(qualityAvg, root=0)
#         hist       = self.comm.bcast(hist, root=0)
#         binEdges   = self.comm.bcast(binEdges, root=0)
            
#         return qualityMin, qualityMax, qualityAvg, hist, binEdges

#     def computeArea(self, *args, **kwargs):
#         mpiPrint("computeArea() has moved to SUmb", comm=self.comm)

#         return 0.0

#     def computeAreaSensitivity(self, *args, **kwargs):
#         mpiPrint("computeAreaSensitivity() has moved to SUmb", comm=self.comm)

#         return []

#     def computeVolume(self, *args, **kwargs):
#         mpiPrint("computeVolume() has moved to SUmb", comm=self.comm)

#         return 0.0

#     def computeVolumeSensitivity(self, *args, **kwargs):
#         mpiPrint("computeVolumeSensitivity() has moved to SUmb",comm=self.comm)

#         return []

#     def releaseMemory(self):
#         """Release all the mesh warping memory"""
#         self.warp.destroyall()
        
#         return 

# # =========================================================================
# #                     Internal Private Functions
# # =========================================================================
        
    def _getFamilyList(self):
        """
        Obtain the family list from fortran. 

        """
        fullList = self.warp.griddata.familylist.transpose().flatten()
        nFamilies = self.warp.griddata.nwallfamilies
        print ('Full list',fullList,nFamilies)
        self.familyList = ["" for i in range(nFamilies)]
        for i in range(nFamilies):
            tempName = fullList[i*32:(i+1)*32]
            newName = ''.join([tempName[j] for j in range(32)]).lower()
            self.familyList[i] = newName.strip()
        # end for

        return

#     def _reOrderIndices(self,  FETopo,  faceBCs,  sym):
#         """This funcion takes the order from self.FE_topo and reorders
#         them according to the Boundary Condition types in each volume
#         class. The sole purpose of this is to facilitate the
#         structural mesh warping algorithim for matrix assembly."""
#         # We want the global indicies ordered according to:
#         #[ freedof ]
#         #[ constrained dof ]

#         ptDOF = np.zeros((FETopo.nGlobal, 3), 'intc')
#         wallDOF = np.zeros((FETopo.nGlobal, 3), 'intc')
#         for ii in range(FETopo.nGlobal):
#             for jj in range(len(FETopo.g_index[ii])):

#                 ivol = FETopo.g_index[ii][jj][0]
#                 i    = FETopo.g_index[ii][jj][1]
#                 j    = FETopo.g_index[ii][jj][2]
#                 k    = FETopo.g_index[ii][jj][3]

#                 N = FETopo.l_index[ivol].shape[0]
#                 M = FETopo.l_index[ivol].shape[1]
#                 L = FETopo.l_index[ivol].shape[2]
              
#                 ptType, number, index1, index2 = \
#                     geo_utils.indexPosition3D(i, j, k, N, M, L)

#                 checkFaces = []
#                 if ptType == 0:
#                     ptDOF[ii] = [0, 0, 0]
#                 else: # Face
#                     if ptType == 1: # Face
#                         checkFaces.append(number)
#                     elif ptType == 2: # Edge
#                         if number in [0, 1, 2, 3]:
#                             checkFaces.append(0)
#                         if number in [4, 5, 6, 7]:
#                             checkFaces.append(1)
#                         if number in [2, 6, 8, 10]:
#                             checkFaces.append(2)
#                         if number in [3, 7, 9, 11]:
#                             checkFaces.append(3)
#                         if number in [0, 4, 8, 9]:
#                             checkFaces.append(4)
#                         if number in [1, 5, 10, 11]:
#                             checkFaces.append(5)
#                     elif ptType == 3: # Corner
#                         if number == 0:
#                             checkFaces.extend([0, 2, 4])
#                         elif number == 1:
#                             checkFaces.extend([0, 3, 4])
#                         elif number == 2:
#                             checkFaces.extend([0, 2, 5])
#                         elif number == 3:
#                             checkFaces.extend([0, 3, 5])
#                         elif number == 4:
#                             checkFaces.extend([1, 2, 4])
#                         elif number == 5:
#                             checkFaces.extend([1, 3, 4])
#                         elif number == 6:
#                             checkFaces.extend([1, 2, 5])
#                         elif number == 7:
#                             checkFaces.extend([1, 3, 5])
#                     # end if
                    
#                     # We now know all faces a point that belong to a
#                     # pt, check each for boundary conditions

#                     for iii in range(len(checkFaces)):
#                         iface = checkFaces[iii]
#                         if faceBCs[ivol][iface] == 1: # Wall
#                             ptDOF[ii] = [1, 1, 1]
#                             wallDOF[ii] = [1, 1, 1]

#                         if faceBCs[ivol][iface] == 2: # Farfield
#                             ptDOF[ii] = [1, 1, 1]
                           
#                         if faceBCs[ivol][iface] == 3:
#                             # Only set it as a symmetry plane if
#                             # nothing is already set
#                             indexSet = np.where(sym==-1)[0][0]
#                             ptDOF[ii][indexSet] = 1

#                         # end if
#                     # end for
#                 # end if
#             # end for
#         # end for
      
#         nus = int(sum(sum(ptDOF)))
#         nuu = FETopo.nGlobal*3-nus
#         nwall = int(sum(sum(wallDOF)))
#         mpiPrint('  -> Total DOF  : %d'%(FETopo.nGlobal*3), comm=self.comm)
#         mpiPrint('  -> Unknown DOF: %d'%(nuu), comm=self.comm)
#         mpiPrint('  -> Known DOF  : %d'%(nus), comm=self.comm)
#         mpiPrint('  -> Wall DOF   : %d'%(nwall), comm=self.comm)
#         # We will forgo the g_index reorganization...it is not
#         # strictly necessary We want l_index[ivol] to be of size
#         # (nu,nv,nw,3) with each entry pointing to the dof in the
#         # global matrix

#         freeDOFCount = 0
#         constrDOFCount = 0
#         l_index = []
#         for ivol in range(len(FETopo.l_index)):
#             l_index.append(np.zeros((3,
#                                   FETopo.l_index[ivol].shape[0],
#                                   FETopo.l_index[ivol].shape[1],
#                                   FETopo.l_index[ivol].shape[2]),
#                                  'intc', order='f'))
#         # end for
            
#         for ii in range(FETopo.nGlobal):
#             for iii in range(3):
#                 if ptDOF[ii][iii] == 0:
#                     for jj in range(len(FETopo.g_index[ii])):
#                         ivol = FETopo.g_index[ii][jj][0]
#                         i    = FETopo.g_index[ii][jj][1]
#                         j    = FETopo.g_index[ii][jj][2]
#                         k    = FETopo.g_index[ii][jj][3]
#                         l_index[ivol][iii, i, j, k] = freeDOFCount
#                     # end for
#                     freeDOFCount += 1
#                 # end if
#                 if ptDOF[ii][iii] == 1:
#                     for jj in range(len(FETopo.g_index[ii])):
#                         ivol = FETopo.g_index[ii][jj][0]
#                         i    = FETopo.g_index[ii][jj][1]
#                         j    = FETopo.g_index[ii][jj][2]
#                         k    = FETopo.g_index[ii][jj][3]
#                         l_index[ivol][iii, i, j, k] = nuu + constrDOFCount
#                     # end for
#                     constrDOFCount += 1
#                 # end if
#             # end for (iii loop)
#         # end for (ii loop)

#         # Lastly, we need to flatten the l_index for fortran use

#         l_index_flat = []
#         lPtr = [0] # -> Zero Based Here
#         lSizes = np.zeros((len(l_index), 3), 'intc')
#         for i in range(len(l_index)):
#             l_index_flat.extend(l_index[i].flatten('f'))
#             lPtr.append(lPtr[-1] + l_index[i].size)
#             lSizes[i] = [l_index[i].shape[1],
#                          l_index[i].shape[2],
#                          l_index[i].shape[3]]
#         # end for

#         return nuu, nus, l_index_flat, lPtr, lSizes

#     def _checkSizes(self, sizes, blockDims):
#         """
#         Make sure the sizes we specify are not larger than the
#         dimension of the block. This can easily happen with "2d"
#         grids, ie. 3d grids only 2 nodes wide

#         We also need to check that odd blockDims have off sizes and
#         even blockDims have even sizes.

#         """
#         for i in range(len(sizes)):
#             for j in range(3):

#                 if sizes[i, j] > blockDims[i, j]:
#                     mpiPrint(' * Warning: Reduced n to %d for an edge'\
#                                  %(blockDims[i, j]))
#                     sizes[i, j] = blockDims[i, j]
#                 # end if

#                 sizeEven = np.mod(sizes[i, j], 2)
#                 dimsEven = np.mod(blockDims[i, j], 2)

#                 if sizeEven == 1 and dimsEven == 1:
#                     pass
#                 elif  sizeEven == 0 and dimsEven == 0:
#                     pass
#                 elif sizes[i, j] == 2:
#                     pass
#                 else:
#                     # Decrease size by 1 to make it even
#                     sizes[i, j] -= 1
#                     sizes[i, j] = np.max([sizes[i, j], 2])
                    
#                     mpiPrint(' * Warning: Reduced n to %d for an edge'\
#                                  %(sizes[i, j]))
#                 # end if
#             # end for
#         # end for

#         return sizes

    def _getIndices(self, groupName):
        """
        Try to see if groupName is already set, if not raise an
        exception

        """
        try: 
            indices = self.familyGroup[groupName]['indices']
        except:
            mpiPrint('+----------------------------------------------------+')
            mpiPrint('Error: %s has not been set using addFamilyGroup.'\
                         %(groupName), comm=self.comm)
            mpiPrint('+----------------------------------------------------+')
        # end try

        return indices

#     def _getCGNSCoords(self):
#         """
#         Internal function to return the block coordinates and
#         Bcoundary conditions. It used in several places hence the the function

#         """

#         nBlocks = self.warp.block.nblocktotal
#         dimsTemp = self.warp.blockdims(nBlocks).reshape((nBlocks, 3))
#         dimsTemp = self.comm.bcast(dimsTemp, root=0)

#         temp0, temp1, CGNSIDs = self.warp.warpmeshsolid_preprocess(
#             nBlocks)

#         # CGNSIDs is 1 based, make it zero based:
#         CGNSIDs -= 1

#         # Reshpae coordinates
#         temp0 = temp0.reshape((nBlocks, 26, 3))
#         temp1 = temp1.reshape((nBlocks, 6))
 
#         # Re-order coordiantes based on CGNSIDs to produce
#         # "standard" original CGNS ordering
#         coords = np.zeros_like(temp0)
#         BCs    = np.zeros_like(temp1)
#         blockDims = np.zeros_like(dimsTemp)
#         for i in range(nBlocks):
#             coords[CGNSIDs[i]] = temp0[i]
#             BCs[CGNSIDs[i]] = temp1[i]
#             blockDims[CGNSIDs[i]] = dimsTemp[i]
#         # end for

#         return nBlocks, coords, BCs, blockDims
      
    def _checkOptions(self, solverOptions):
        """
        Check the solver options against the default ones
        and add opt
        ion iff it is NOT in solverOptions

        """
        for key in self.solverOptionsDefault.keys():
            if not(key in solverOptions.keys()):
                solverOptions[key] = self.solverOptionsDefault[key]
            else:
                self.solverOptionsDefault[key] = solverOptions[key]	
            # end if
        # end for
        #print('solveroptions',solverOptions)
        return solverOptions
