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
import pywarpustruct
from petsc4py import PETSc

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
    def __init__(self, optionsList, bgFileList=[], comm=None):
        """
        Create the MultiUSMesh object.

        INPUTS:

        optionsList: list of dictionaries -> List containing dictionaries that will
        be used to initialize multiple pyWarp instances. Each dictionary must have
        a structured CGNS filename in its 'gridFile' field.
        ATTENTION: The CGNS files should follow the same ordering as the combined
        CGNS file given to ADflow.

        bgFileList: list of strings -> List containing file names of CGNS files that
        contains just background meshes, and were combined to make the full CGNS file
        given to ADflow
        ATTENTION: The CGNS files should follow the same ordering as the combined
        CGNS file given to ADflow.

        Ney Secco 2017-02
        """

        # Assign communicator if we do not have one yet
        if comm is None:
            self.comm = MPI.COMM_WORLD
        else:
            self.comm = comm

        # Get rank of the current proc
        myID = self.comm.Get_rank()

        #------------------------------------------------------
        # READING NEARFIELD MESHES (The ones that will be warped)
        #

        # Initialize list of pyWarp instances
        self.meshes = []

        # Initialize list that will contain indices to slice the global surface
        # coordinates vector into the pieces that belong to each instance
        self.surfSliceIndices = [0]

        # Loop over the multiple CGNS files to initialize the corresponding pyWarp instances.
        # We use the same loop to count how many surface nodes we have in each instance.
        for meshOptions in optionsList:

            # Initialize a pyWarp instance with the current options
            currMesh = pywarpustruct.USMesh(options=meshOptions, comm=self.comm)

            # Append the instance to the list
            self.meshes.append(currMesh)

            # Get the points in the current pyWarp instance.
            # Note that this will force pyWarp to build its own connectivities, which
            # will (luckily) be the same structure seen by ADflow.
            pts = currMesh.getSurfaceCoordinates()

            # Get the number of nodes that belong to the current instance
            numSurfPts = pts.shape[0]

            # Add index to help slice the global surface array later on
            self.surfSliceIndices.append(self.surfSliceIndices[-1] + numSurfPts)

        # Get list of nodes in the current proc, divided by instance.
        # That is, volNodesList[i] has the volume nodes in the current proc that
        # belong to the i-th pyWarp instance.
        # The same function also returns the total number of coordinates (numCoor = 3*numNodes)
        # of all volume meshes.
        volNodesList, numCoorTotal = self.getWarpGrid()

        # Store the total number of coordinates
        self.numCoorTotal = numCoorTotal

        #------------------------------------------------------
        # GENERATING pyWarpMulti GLOBAL VOLUME VECTOR FOR NEARFIELD MESHES
        #
        # We loaded all meshes nearfield meshes (the ones that will be warped)
        # Now it is time to initialize the global PETSc vector that will hold
        # all nearfield volume nodes following the CGNS ordering.
        #
        # In the end, we need a PETSc global vector with the following structure:
        #
        # |------p0--------|----------p1---------|------p2-------|
        # xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx = warpVolPtsVec
        # |------mesh 0------|------mesh 1------|----mesh 2------|
        # |--p0--|-p1-|--p2--|--p0--|--p1--|-p2-|-p0-|-p1-|--p2--| p:proc
        #
        # The upper part shows where the elements are stored, while the
        # lower shows where the elements are coming from.
        #
        # volNodesList is a local list that contains the volume nodes in this proc,
        # divided per instance.
        # For example, if I am in proc 1 and call volNodesList[0], I get all
        # elements covered by the "p1" tag, under the "mesh 0 tag".
        #
        # Note that instead of number of nodes, we work with number of degrees of freedom, so
        # that we can represent volume coordinates as a flattened array.
        # Remember that numDOF = numVolNodes*3
        #
        # Fortunately, pyWarp does not duplicate volume nodes, so we can use the number of nodes
        # itself to determine which slice of the CGNS file that belongs to the current proc.
        # This means that the automatic indexList generated by fillGlobalVector will do the trick!

        # Are you ready for some PETSc?

        # Initialize a PETSc vector
        warpVolPtsVec = PETSc.Vec()

        # Here we set the global vector size, and let PETSc do the partitioning.
        # Fortunately, the variable instanceOffset has the total number of coordinate points.
        warpVolPtsVec.createMPI(size=numCoorTotal)

        # Call the function that will fill all entries of the global vector.
        # We will store the automatically generated slicing indices for the future
        # volume mesh updates.
        self.volIndexList = fillGlobalVector(warpVolPtsVec, volNodesList, comm)

        # Just for the record, the i-th entry of the volIndexList is a 1x2 array that contains
        # the first and last global vector indices that correspond to the i-th pyWarp instance
        # coordinates that are in the current proc.

        # Store the global vector with the nearfield mesh volume nodes.
        self.warpVolPtsVec = warpVolPtsVec

        #------------------------------------------------------
        # READING BACKGROUND MESHES
        #
        # We will initialize a global vector that will store all background mesh nodes.
        # It will be analogous to the nearfield global vector.
        #
        # First we will load each background mesh, to get the volume nodes,
        # then we can initialize the PETSc vector and assign the corresponding slices.

        # Initialize list to hold volume nodes (in the current proc) of all instances.
        # That is, volNodesList[i] gives the volume nodes of the i-th pyWarp instance that belong
        # to the current proc.
        volNodesList = []

        # Initialize counter to store the total number of coordinates (numCoor = 3*numNodes)
        # of all volume mesh.
        numCoorTotal = 0

        for bgFile in bgFileList:

            #=========================================================#
            # THIS IS A MESSY (HOPEFULLY TEMPORARY) WAY OF LOADING THE
            # BACKGROUND MESH NODES. IF YOU COME UP WITH A BETTER WAY
            # TO GET volNodes, PLEASE ADD IT HERE.
            # volNodes is a flattened vector that contains the background
            # mesh volume nodes that belong to the current proc.

            # Let's try using pyWarp's CGNS loader to extract the bakground nodes.
            # However, we will have to trick pyWarp in order to do this, since it
            # expects a surface mesh in the file.
            # So we will make a copy of the background mesh file, assign an arbitrary
            # wall surface, and then load it with pyWarp

            # Only the root proc will modify the input file
            if myID == 0:

                # Make a copy of the background mesh file
                os.system('cp '+bgFile+' tmp_bg_file.cgns')

                # Create a temporary BC file
                with open('tmp_bcdata.dat','w') as fid:
                    fid.write('1 iLow BCwall wall\n')

                # Use CGNS utils to modify the BCs
                os.system('cgns_utils overwritebc tmp_bg_file.cgns tmp_bcdata.dat')

            # Create dummy set of options just to load the CGNS file
            optionsDict = {
                'gridFile':'tmp_bg_file.cgns',
                'warpType':'unstructured',
            }

            # Initialize a pyWarp instance with the current options
            currMesh = pywarpustruct.USMesh(options=optionsDict, comm=self.comm)

            # The root proc can remove the temporary files
            if myID == 0:

                # Make a copy of the background mesh file
                os.system('rm tmp_bg_file.cgns')
                os.system('rm tmp_bcdata.dat')

            # Dummy call to initialize warping
            currMesh._setInternalSurface()

            # Get volume nodes.
            # volNodes is a flattened vector that contains the background
            # mesh volume nodes that belong to the current proc.
            volNodes = currMesh.getWarpGrid()

            #=========================================================#

            # Store the nodes of the current instance in the list
            volNodesList.append(volNodes)

            # Get number of coordinates on the current processor, for the current pyWarp instance.
            numCoor = len(volNodes)

            # Each processor needs to know how many coordinates are in the other processors.
            # We use an allreduce operation to sum the number of coordinates in all procs
            numCoor_all = comm.allreduce(numCoor, MPI.SUM)

            # Increment counter for the next instance
            numCoorTotal = numCoorTotal + numCoor_all

        #------------------------------------------------------
        # GENERATING pyWarpMulti GLOBAL VOLUME VECTOR FOR THE BACKGROUND MESHES

        # Whew... We loaded all meshes...
        # Now it is time to initialize the global PETSc vector.
        #
        # In the end, we need a PETSc global vector with the following structure:
        #
        #
        # |------p0--------|----------p1---------|------p2-------|
        # xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx = warpVolPtsVec_BG
        # |------mesh 0------|------mesh 1------|----mesh 2------|
        # |--p0--|-p1-|--p2--|--p0--|--p1--|-p2-|-p0-|-p1-|--p2--| p:proc
        #
        # The upper part shows where the elements are stored, while the
        # lower shows where the elements are coming from.
        #
        # Here is a very important thing to keep in mind. Even though a node was read
        # in proc i, it may end up in the global vector slice that belongs to proc j.
        # We allow this because it is important to keep the pyWarp INSTANCES contiguous
        # in the global vector.

        # Here comes PETSc again...

        # Initialize a PETSc vector
        warpVolPtsVec_BG = PETSc.Vec()

        # Here we set the global vector size, and let PETSc do the partitioning.
        # Fortunately, the variable instanceOffset has the total number of coordinate points.
        warpVolPtsVec_BG.createMPI(size=numCoorTotal)

        # Now we will populate the vector with the nodal coordinates.
        # We do not need to store the slicing indices here since the background mesh
        # will not be warped.
        fillGlobalVector(warpVolPtsVec_BG, volNodesList, comm)

        # Store the global vector with the background mesh volume nodes.
        self.warpVolPtsVec_BG = warpVolPtsVec_BG

        #------------------------------------------------------
        # Initialize fields for other PETSc objects that will be created later on
        self.surfScatter = None
        self.volScatter = None

    def getSurfaceCoordinates(self):
        """Returns all defined surface coordinates on this processor, with the ADflow ordering

        REMEMBER THAT PETSc CREATES ITS VECTORS BY REFERENCE, NOT COPY!!!
        So if you do:

        a = self.getSurfaceCoordinates()
        self.setSurfaceCoordinates(b)
        
        Then you will get a==b in the end!

        Returns
        -------
        pts : numpy array size (N,3)
            Specified surface coordinates residing on this
            processor. This may be empty array, size (0,3)

        Ney Secco 2017-02
        """

        # Here we need to do a reverse scatter to take values from the pyWarp surface vector
        # an pass them back to the ADflow surface vector format

        # Do the reverse scatter (mode=True) if it is defined
        if self.surfScatter is None:
            if self.comm.Get_rank() == 0:
                raise NameError('The pyWarp-ADflow scatter is not initialized. Run self.setSurfaceDefinition first.')
        elif self.surfScatter != 'noScatter':
            self.surfScatter.scatter(self.warpSurfPtsVec, self.surfPtsVec, addv=None, mode=True)

        # Now get the slice of the global ADflow surface vector that belongs to the current proc
        # Remember to reshape it since PETSc is working with vectors
        pts = self.surfPtsVec.getArray().reshape((-1,3))

        return pts

    def getSurfaceCoordinates_pyWarp(self):
        """Returns all defined surface coordinates on this processor, with the pyWarp ordering

        Returns
        -------
        pts : numpy array size (N,3)
            Specified surface coordinates residing on this
            processor. This may be empty array, size (0,3)

        Ney Secco 2017-02
        """

        # Initialize array to hold all surface points
        pts = np.zeros((self.surfSliceIndices[-1],3))

        # Initialize instance counter
        instanceID = 0

        # Loop over every mesh object to get a slice of the surface point array
        for mesh in self.meshes:

            # Get coordinates of the current instance
            curr_pts = mesh.getSurfaceCoordinates()
            
            # Assign current set of coordinates to the correct slice of the global array
            pts[self.surfSliceIndices[instanceID]:self.surfSliceIndices[instanceID+1],:] = curr_pts

            # Increment instance counter
            instanceID = instanceID + 1

        return pts

    def setSurfaceCoordinates(self, pts):
        """Sets all surface coordinates on this processor, with pts given in ADflow ordering

        Parameters
        ----------
        pts : numpy array, size(N, 3)
            The coordinate to set. This MUST be exactly the same size as
            the array obtained from getSurfaceCoordinates()

        Ney Secco 2017-02
        """

        # We need to take these user-provided points and send them to the PETSc global vector.
        # We do this because the PETSc vector will be used in the scatter procedure.
        self.surfPtsVec.setArray(pts)

        # Here we need to do a forward scatter to take values from the ADflow surface vector
        # an pass them back to the pyWarp surface vector format

        # Do the forward scatter (mode=None) if it is defined
        if self.surfScatter is None:
            if self.comm.Get_rank() == 0:
                raise NameError('The pyWarp-ADflow scatter is not initialized. Run self.setSurfaceDefinition first.')
        if self.surfScatter != 'noScatter':
            self.surfScatter.scatter(self.surfPtsVec, self.warpSurfPtsVec, addv=None, mode=None)

        # Now we can get the pyWarp surface coordinates tha belong to the current proc
        # Remember to reshape it since PETSc is working with vectors
        warpPts = self.warpSurfPtsVec.getArray().reshape((-1,3))

        # Then we can use the equivalent function that uses pyWarp ordering
        self.setSurfaceCoordinates_pyWarp(warpPts)

    def setSurfaceCoordinates_pyWarp(self, pts):
        """Sets all surface coordinates on this processor, with pts given in pyWarp ordering

        Parameters
        ----------
        pts : numpy array, size(N, 3)
            The coordinate to set. This MUST be exactly the same size as
            the array obtained from getSurfaceCoordinates_pyWarp()

        Ney Secco 2017-02
        """

        # Initialize instance counter
        instanceID = 0

        # Loop over every mesh object to get a slice of the surface point array
        for mesh in self.meshes:

            # Get set of points that belongs to the current instance
            curr_pts = pts[self.surfSliceIndices[instanceID]:self.surfSliceIndices[instanceID+1],:]

            # Set the current points to the current instance
            mesh.setSurfaceCoordinates(curr_pts)

            # Increment instance counter
            instanceID = instanceID + 1

    def setExternalMeshIndices(self, ind):
        """
        Set the indices defining the transformation of an external
        solver grid to the original CGNS grid. This is required to use
        USMesh functions that involve the word "Solver" and warpDeriv.
        The indices must be zero-based.

        Parameters
        ----------
        ind : numpy integer array[3*n] where n is the number of nodes stored in the proc.
              The list of indicies this processor needs from the common mesh file.
              ind[3*i:3*i+3] represents the position of each coordinate of the i-th ADflow node
              in the CGNS file.
              ind has 0-based indices.

        Ney Secco 2017-02
        """

        # Here we will do the mapping between pyWarpMulti and ADflow volume nodes.
        # ASSUMPTION: We assume that the pyWarp instances are loaded in the same order
        # as the full CGNS file loaded by ADflow. In other words, if you concatenate
        # all blocks in the pyWarp CGNS files, you will get the same ordering as the
        # combined CGNS file given to ADflow.

        # Here is the main logic of the function:
        # - Imagine that we create a big [N by 3] array with all nodal coordinates in
        #   the CGNS file. Then we flatten it to get a [N*3] array. This global array
        #   is split across multiple procs by pyWarp. Therefore, pyWarp has this global
        #   vector.
        # - ADflow gives us the indices of this global array that the current proc needs (this is ind).
        #   Thus we can use this information to create a PETSc VecScatter operation, since ind is
        #   the PETSc index set itself!
        
        #----------------------------------------------------
        # INITIALIZE PETSc VECTOR

        # Here we initialize a global PETSc vector that will store the volume coordinates in
        # ADflow ordering.

        # Initialize vector
        volPtsVec = PETSc.Vec()

        # For now let's use the indices themselves to set the vector, since ind has the size of
        # the local number of coordinates. THe data itself will be overwritten later on, we just
        # want the corret size for now.
        # When we use Vec.createWithArray, PETSc will concatenate all local vectors
        # to create the global one.
        # We also use np.array to make a copy of ind.
        volPtsVec.createWithArray(np.array(ind,dtype=float))

        # Store this global vector
        self.volPtsVec = volPtsVec

        #----------------------------------------------------
        # NEARFIELD MESH MAPPING

        # Our main goal in this section is to build the Scatter context to send data from
        # ADflow volume nodes to pyWarpMulti volume nodes.
        # pyWarpMulti will not warp the background mesh. Therefore, we will remove from the
        # nearfield scatter to avoid unnecessary communication during the optimization.

        # Find indices of ind that will point to nearfield nodes. We can find them using
        # the number of nearfield coordinates that we stored during the initialization
        adflowIndex = np.where(ind < self.numCoorTotal)[0]

        # Now take the corresponding links. Note that these will point to the CGNS ordering,
        # which is the one used by pyWarp
        warpIndex = ind[adflowIndex]

        # Get the range of the global ADflow volume vector owned by the current proc
        (indexStartADflowPts, indexEndADflowPts) = volPtsVec.getOwnershipRange()

        # Now we need to offset the indices given by the np.where, since they are originally
        # defined in local coordinates. (Note that warpIndex already got information in
        # the CGNS ordering, and does not need to be offset).
        adflowIndex = adflowIndex + indexStartADflowPts

        # Make sure indices are defined as integers
        adflowIndex = np.array(adflowIndex, dtype='int32')

        ### We need to build the index sets that will define the VecScatter operation.

        # This index set indicates which indices will be taken for the ADflow volume vector
        volPtsVecIS = PETSc.IS()
        volPtsVecIS.createGeneral(adflowIndex)

        # This index set indicates which indices will be taken for the pyWarp volume vector
        warpVolPtsVecIS = PETSc.IS()
        warpVolPtsVecIS.createGeneral(warpIndex)

        ### It is time to create the scatter context

        # Initialize VecScatter
        volScatter = PETSc.Scatter()

        # Create Scatter using the index sets
        # We follow this syntax to create the scatter context:
        # Scatter.create(Vec_from, IS_from, Vec_to, IS_to)
        # IS_from and IS_to should have the same size.
        volScatter.create(self.volPtsVec, volPtsVecIS, self.warpVolPtsVec, warpVolPtsVecIS)

        # Store this scatter
        self.volScatter = volScatter

        #----------------------------------------------------
        # BACKGROUND MESH ASSIGNMENT

        # Our main goal in this section is to assign the background volume nodes to the global
        # volume vector in ADflow ordering. We have to do this just once here, since the
        # background meshes will not be updated during the optimization.

        # Find indices of ind that will point to background nodes. We can find them using
        # the number of nearfield coordinates that we stored during the initialization
        adflowIndex = np.where(ind >= self.numCoorTotal)[0]

        # Now take the corresponding links. Note that these will point to the CGNS ordering,
        # which is the one used by pyWarp
        warpIndex = ind[adflowIndex]

        # Since we will take values from a global vector that only has background nodes,
        # we need to offset the pyWarp indices to take into acount that the background nodes
        # are after all nearfield nodes
        warpIndex = warpIndex - self.numCoorTotal

        # Now we need to offset the indices given by the np.where, since they are originally
        # defined in local coordinates. (Note that warpIndex already got information in
        # the CGNS ordering, and does not need to be offset).
        adflowIndex = adflowIndex + indexStartADflowPts  

        # Make sure indices are defined as integers
        adflowIndex = np.array(adflowIndex, dtype='int32')

        ### We need to build the index sets that will define the VecScatter operation.

        # This index set indicates which indices will be taken for the ADflow volume vector
        volPtsVecIS = PETSc.IS()
        volPtsVecIS.createGeneral(adflowIndex)

        # This index set indicates which indices will be taken for the pyWarp volume vector
        warpVolPtsVecIS = PETSc.IS()
        warpVolPtsVecIS.createGeneral(warpIndex)

        ### It is time to create the scatter context

        # Initialize VecScatter
        volScatter = PETSc.Scatter()

        # Create Scatter using the index sets
        # We follow this syntax to create the scatter context:
        # Scatter.create(Vec_from, IS_from, Vec_to, IS_to)
        # IS_from and IS_to should have the same size.
        volScatter.create(self.volPtsVec, volPtsVecIS, self.warpVolPtsVec_BG, warpVolPtsVecIS)

        # Now let's use this scatter context right away to send the background nodes from
        # the pyWarp-ordered global vector the the ADflow-ordered global vector.
        # In this case we use the reverse scatter (mode=True).
        volScatter.scatter(self.warpVolPtsVec_BG, self.volPtsVec, addv=None, mode=True)

        ### DESTRUCTION

        # Now we can erase most PETSc objects and variables associated with the background mesh
        # since they are no longer needed during the optimization
        self.warpVolPtsVec_BG.destroy()
        self.warpVolPtsVec_BG = 'used'
        volScatter.destroy()
        volPtsVecIS.destroy()
        warpVolPtsVecIS.destroy()

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

        # Hopefully, the CGNS ordering will save us.
        # The automatic connectivities generated within pyWarp might match
        # ADflow if we load the separate CGNS files in order as well!

        # First we get the updated volume nodes from all instances
        # Remember that volNodesList[i] has the volume nodes in the current proc that
        # belong to the i-th pyWarp instance.
        volNodesList, numCoorTotal = self.getWarpGrid()

        # Now we use the list to update the global PETSc vector, in pyWarp ordering, for
        # the scatter operation
        # Note that we can use the slicing intervals automatically generated during the
        # initialization method (self.volIndexList)
        fillGlobalVector(self.warpVolPtsVec, volNodesList, self.comm, indexList=self.volIndexList)

        # Now let's use this scatter context right away to send the nearfield nodes from
        # the pyWarp-ordered global vector the the ADflow-ordered global vector.
        # In this case we use the reverse scatter (mode=True).
        # We also check if the scatter was initialized
        if self.volScatter is None:
            if self.comm.Get_rank() == 0:
                raise NameError('The pyWarp-ADflow scatter is not initialized. Run self.setExternalMeshIndices first.')
        self.volScatter.scatter(self.warpVolPtsVec, self.volPtsVec, addv=None, mode=True)

        # Get the local coordinates from the PETSc global vector
        solverGrid = self.volPtsVec.getArray()

        return solverGrid

    def setSurfaceDefinition(self, pts, conn=None, faceSizes=None, distTol=1e-6):
        """This is the master function that determines the definition of the
        surface to be used for the mesh movement. This surface may be
        supplied from an external solver (such as ADflow) or it may be
        generated by pyWarpUStruct internally.

        The main goal of this function is the creation of the PETSc objects
        to handle the data passing between ADflow and pyWarp. These objects
        will be stored in self.ptsVec (ADflow global surface coordinate vector),
        self.warpPtsVec (pyWarp global surface coordinate vector), and self.surfScatter
        (scatter context to communicate both vectors). A "global" vector means that it
        is distributed across many procs.

        I also want to outline that there is a big difference between the global surface
        vectors and the global volume vectors. The global volume vector preserves the
        original CGNS ordering (it is instance contiguous), while we rearrange the
        surface vector elements to gather data by proc (it is proc contiguous).
        Here is a better illustration:

        global volume vector:

        |------p0--------|----------p1---------|------p2-------|
        xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        |------mesh 0------|------mesh 1------|----mesh 2------|
        |--p0--|-p1-|--p2--|--p0--|--p1--|-p2-|-p0-|-p1-|--p2--| p:proc

        global surface vector:

        |------p0--------|----------p1---------|-------p2-------|
        xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        |-m0-|-m1-|--m2--|--m0--|---m1--|--m2--|-m0-|-m1-|--m2--| m:mesh
        |-------p0-------|---------p1----------|-------p2-------|

        The upper part of each drawing shows where the elements are stored, while the
        lower shows where the elements are coming from.
        
        We can rearrange the surface nodes because we do a distance-based mapping, while
        the volume nodes are mapped based on the original CGNS ordering.

        Parameters
        ----------
        pts : array, size (M, 3)
            Nodes on this processor
        conn : int array, size  (sum(faceSizes))
            Connectivity of the nodes on this processor
            THIS IS NOT USED IN THIS FUNCTION, BUT WE KEEP IT TO HAVE THE
            SAME API AS A SINGLE PYWARP INSTANCE
        faceSizes : int Array size (N)
            Treat the conn array as a flat list with the faceSizes giving
            connectivity offset for each element.
            THIS IS NOT USED IN THIS FUNCTION, BUT WE KEEP IT TO HAVE THE
            SAME API AS A SINGLE PYWARP INSTANCE
        distTol: Distance tolerance to flag that a given surface node does not
                 belong to the current pyWarp surface definition in the current proc.

        Ney Secco 2017-02
        """

        # IMPORTS
        from scipy.spatial import cKDTree

        # We will throw the connectivities in the trash!!!
        # We will receive all surface nodes coming from ADflow and map them to the surface
        # nodes we have in pyWarp.
        #
        # Here is the main logic of this function
        # - Each proc receives a set of surface points from ADflow
        # - We do an allGather to distribute all ADflow surface nodes to all procs.
        #   Therefore, each proc will have the global surface vector.
        # - Each proc will take its pyWarp surface nodes and build a KDTree
        # - We use the KDtree to map the pyWarp surface nodes currently owned by this proc
        #   to the indices of the ADflow global surface vector.
        # - We use this mapping to create a petsc scatter context. Then we can use this same context
        #   in the next iterations to handle all communications directly!

        # Create the ADflow global surface vector on every proc
        pts_all = np.vstack(self.comm.allgather(pts))

        # Get the surface nodes currently defined in our pyWarp instances
        warpPts = self.getSurfaceCoordinates_pyWarp()

        # Get the number of pyWarp surface nodes
        numWarpPts = warpPts.shape[0]

        #---------------------------
        ### ADflow-pyWarp Mapping

        # Here we use a KDTree to match the ADflow nodes with the pyWarp nodes
        
        # We need to be careful because some nodes in the pyWarp vector may be repeated (shared by multiple blocks),
        # and we can't leave these repeated nodes out of the mapping. So we allow the KDTree to search for multiple
        # candidates, then we can assign all possible repeated mappings.
        # Here we set the maximum number we expect that a node may be repeated.
        maxRep = 6

        # If the proc has no pyWarp surface nodes, we don't need to create any
        # mapping.
        if numWarpPts != 0:

            # Each proc creates its KDtree using its pyWarp points
            tree = cKDTree(warpPts)

            # Now use the KDTree to find which index from warpPts is correlated to a given
            # node from pts_all.
            # That is pts_all[ii] = warpPts[indexMap[ii]]
            # If a given ADflow node does not match with any node in the tree, then its
            # indexMap value will be numWarpPts (which is out of the range of warpPts)
            # We also allow the KDTree to search for the best k candidates, so we can
            # take care of repeated nodes.
            dist, indexMap = tree.query(pts_all, distance_upper_bound=distTol, k=maxRep)

            # Convert indexMap to a numpy array to facilitate next operations
            indexMap = np.array(indexMap)

            # At this point, indexMap is [numPtsAll x maxRep]. Therefore,
            # indexMap[i,j] gives the j-th candidate index in warpPts that corresponds
            # to the i-th node in pts_all.
            # So now we need to analyze column by column to take into account the multiple
            # candidates.

            # First let's initialize 1D arrays that we will increment for every candidate
            # analysis. Please see the comments over indexPtsCurr and indexWarpPtsCurr to
            # understand the role of these arrays
            indexPts = []
            indexWarpPts = []

            for candID in range(maxRep):

                # Now we need to remove entries from the indexMap that correspond to ADflow nodes
                # that did not match with any pyWarp node in this proc. We can find these nodes
                # because indexMap[ii,candID]==numWarpPts.
                # Let's say we have numWarpPts = 5, and our current indexMap is
                # indexMap[:,candID] = [0 1 5 3 2 5 4]
                # This means that the 0th pts node is linked to the 0th warpPts node,
                # the 1st pts node is linked to the 1st warpPts node, the 2nd pts node
                # did not find any match in warpPts, the 3rd pts node is linked to the 3rd warpPts
                # node, the 4th pts node is linked to the 2nd warpPts node, and so on...
                # In the end, the important index relations are:
                #  pts -> warpPts
                #  0th -> 0th
                #  1th -> 1th
                #  3rd -> 3rd
                #  4th -> 2nd
                #  6th -> 4th
                # So we will create two arrays to express these relationships.
                # The first column will be indexPtsCurr, and the second one will be indexWarpPtsCurr
                # These arrays will be concatenated into indexPts and indexWarpPts to gather results
                # for all candidate orders.

                # Find indices of indexMap that are smaller than numWarpPts. These are
                # the indices that managed to find a match in the current proc.
                indexPtsCurr = np.where(indexMap[:,candID] < numWarpPts)[0]

                # Now take the corresponding links
                indexWarpPtsCurr = indexMap[indexPtsCurr,candID]

                # Concatenate these arrays into the major arrays that gather results for all
                # candidates
                indexPts = np.hstack([indexPts,indexPtsCurr])
                indexWarpPts = np.hstack([indexWarpPts,indexWarpPtsCurr])

            # We don't need the local copy of the global surface vector anymore.
            # We also don't need the unsorted index mapping as well.
            del pts_all
            del indexMap

        else:

            # Just create an empty mapping
            indexPts = np.array([])
            indexWarpPts = np.array([])

        # Are you ready?
        # The fun PETSc stuff starts here...

        #---------------------------
        # Now we will create a global PETSc vector with all the ADflow surface coordinates
        # from each proc
        ptsVec = PETSc.Vec()

        # Each proc will contribute with its surface nodes to build the global vector
        # Note that we will work with a flattened surface vector
        ptsVec.createWithArray(pts.flatten())

        # Now we need to create another PETSc vector to hold the pyWarp surface coordinates
        warpPtsVec = PETSc.Vec()

        # Each proc will contribute with its surface nodes to build the global vector
        # Note that we will work with a flattened surface vector
        warpPtsVec.createWithArray(warpPts.flatten())

        #---------------------------
        # We need to get ready for the VecScatter operation, in which we will take values
        # from the ADflow surface coordinates and set them into pyWarp coordinates.
        # Just to remember, a VecScatter operation will take values from a parallel vectors
        # and place then into another parallel vector, which may have different partitionings.

        # Get the range of the global pyWarp surface vector owned by the current proc
        (indexStartWarpPts, indexEndWarpPts) = warpPtsVec.getOwnershipRange()

        # Now we need to offset the indices given by the KDTree, since they are originally
        # defined in local coordinates. (Note that indexPts was built already in global coordinates
        # since we used the global vector in the KDTree).
        indexWarpPts = indexWarpPts + indexStartWarpPts

        # Up to now, all mapping indices were defined for the [Nx3] arrays, but the PETSc vectors are [3*N x 1],
        # so we need to convert them.
        indexPts1D = convertMapping2Dto1D(indexPts,3)
        indexWarpPts1D = convertMapping2Dto1D(indexWarpPts,3)

        #---------------------------
        ### We need to build the index vectors that will define the VecScatter operation.

        # This index set indicates which indices will be taken for the ADflow surface vector
        ptsVecIS = PETSc.IS()
        ptsVecIS.createGeneral(indexPts1D)

        # This index set indicates which indices will be taken for the pyWarp surface vector
        warpPtsVecIS = PETSc.IS()
        warpPtsVecIS.createGeneral(indexWarpPts1D)

        #---------------------------
        ### It is time to create the scatter context

        # Initialize VecScatter
        surfScatter = PETSc.Scatter()

        # Create Scatter using the index sets
        # We follow this syntax to create the scatter context:
        # Scatter.create(Vec_from, IS_from, Vec_to, IS_to)
        # IS_from and IS_to should have the same size.
        surfScatter.create(ptsVec, ptsVecIS, warpPtsVec, warpPtsVecIS)

        #---------------------------
        # Finally we can store the PETSc objects that will be needed later on to match surface information
        self.surfPtsVec = ptsVec
        self.warpSurfPtsVec = warpPtsVec
        self.surfScatter = surfScatter

    def getWarpGrid(self):
        """
        Return the current grids. This function is typically unused. See
        getSolverGrid for the more useful interface functionality.

        This only returns the nearfield meshes.

        Returns
        -------
        volNodesList, list of 1D numpy arrays, real: These are the local volume nodes (in a flat 1D array)
        of each instance. That is, volNodesList[i] has the volume nodes stored in the local proc for
        the i-th pyWarp instance.

        numCoorTotal: The total number of coordinates, across all procs and pyWarp instances.

        Ney Secco 2017-02
        """

        # Initialize list to hold volume nodes (in the current proc) of all instances.
        # That is, volNodesList[i] gives the volume nodes of the i-th pyWarp instance that belong
        # to the current proc.
        volNodesList = []

        # Initialize counter to store the total number of coordinates (numCoor = 3*numNodes)
        # of all volume mesh.
        numCoorTotal = 0

        # Loop over the multiple CGNS files to initialize the corresponding pyWarp instances
        for currMesh in self.meshes:

            # Get volume nodes.
            # volNodes is a flattened vector that contains the background
            # mesh volume nodes that belong to the current proc.
            volNodes = currMesh.getWarpGrid()

            # Store the nodes of the current instance in the list
            volNodesList.append(volNodes)

            # Get number of coordinates on the current processor, for the current pyWarp instance.
            numCoor = len(volNodes)

            # Each processor needs to know how many coordinates are in the other processors.
            # We use an allreduce operation to sum the number of coordinates in all procs
            numCoor_all = self.comm.allreduce(numCoor, MPI.SUM)

            # Increment counter for the next instance
            numCoorTotal = numCoorTotal + numCoor_all

        # RETURNS
        return volNodesList, numCoorTotal

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
        This calls the mesh warping method for each pyWarp instance.

        This will update the volume coordinates internally in each instance.

        Ney Secco 2017-02
        """

        # Get proc ID
        myID = self.comm.Get_rank()

        # Print log
        if myID == 0:
            print('')
            print('Starting pyWarpMulti mesh warping')

        # Set mesh counter
        meshCounter = 1

        # Loop over all instances
        for mesh in self.meshes:

            # Print log
            if myID == 0:
                print('')
                print(' warping mesh',meshCounter,'of',len(self.meshes))

            # Warp current instance
            mesh.warpMesh()

            # Print log
            if myID == 0:
                print('')
                print(' Done')

            # Increment counter
            meshCounter = meshCounter + 1

        # Print log
        if myID == 0:
            print('')
            print('pyWarpMulti successfully warped all instances!')

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

# =========================================================================
# Other simple functions
# =========================================================================

def convertMapping2Dto1D(mapping2D,numDim):

    '''
    This function takes an index mapping for 2D coordinate vector
    and transforms it into a mapping for an equivalent 1D flattened
    coordinate vector.

    INPUTS:

    mapping2D: array[n] -> Array that has an index mapping for a 2D
    coordinate array

    numDim: integer -> Dimension of the 2D matrix we are mapping to

    OUTPUTS:

    mapping1D: array[3*n] -> Array that has an index mapping for an
    equivalent flattened coordinate array

    Ney Secco 2017-02
    '''

    # Get number of elements and dimension of the matrix
    numElem = len(mapping2D)
    
    # Initialize new mapping array
    mapping1D = np.zeros(numElem*numDim,dtype='int32')

    # Populate indices corresponding to each dimension
    for iDim in range(numDim):
        mapping1D[iDim::numDim] = mapping2D[:]*numDim + iDim

    # Return flattened mapping
    return mapping1D

def fillGlobalVector(globalVec, localList, comm, indexList=None):

    '''
    This function will take data from localList, which contains local vector with
    coordinates from all instances, and send it to the global vector globalVec.

    INPUTS:
    globalVec: PETSc global vector -> Global vector that contains coordinates of all
    instances and all procs. It follows this structure:
    
     |------p0--------|----------p1---------|------p2-------|
     xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx = globalVec
     |------mesh 0------|------mesh 1------|----mesh 2------|
     |--p0--|-p1-|--p2--|--p0--|--p1--|-p2-|-p0-|-p1-|--p2--| p:proc
    
    The upper part shows where the elements are stored, while the
    lower shows where the elements are coming from.

    Here is a very important thing to keep in mind. Even though a node was read
    in proc i, it may end up in the global vector slice that belongs to proc j.
    We allow this because it is important to keep the pyWarp INSTANCES contiguous
    in the global vector.

    localList: list of local 1D arrays -> List containing coordinates that belong to proc 0.
    localList[i] should have a 1D array with the coordinates of pyWarp instance i that
    belong to the current proc. Each element of localList corresponds to a lower "p[i]" tag
    on the ASCII art.

    indexList: list of list[2] -> List containing the initial and final indices of the global
    vector that define where the local data should be placed. For instance, if indexList[i]=[a,b],
    then we will take all data in localList[i] and place it into globalVec[a:b].
    If the user provide None, this function will generate the ordering on its own, assuming the
    structure shown in the drawing. Then the user can store the returned indexList for future use.

    OUTPUTS:

    indexList: list of list[2] -> This is will be list provided by the user, or the one that was
    automatically generated by the function.

    This function also modifies globalVec.

    Ney Secco 2017-02
    '''

    # Get current proc ID
    myID = comm.Get_rank()

    # Generate an indexList if the user provided None
    if indexList is None:

        # Initialize list that will contain indices to slice the global vector.
        # The i-th entry of the list is a 1x2 array that contains the first and last
        # indices that belong to the current proc in the i-th pyWarp instance
        indexList = []

        # Initialize offset variable to take into account how many elements you have in
        # each instance
        instanceOffset = 0

        for instanceID in range(len(localList)):

            # Get number of degrees of freedom on the current processor, for the current pyWarp instance.
            numDOF = len(localList[instanceID])

            # Each processor needs to know how many nodes are in the other processors.
            # We use an allGather operation to send and receive all number of DOF.
            numDOF_all = np.hstack(comm.allgather([numDOF]))

            # Create array to store the interval that belongs to the current proc
            # at the current pyWarp instance. myIndexList will store the first and
            # last index of the current slice
            myIndexList = [0,0]

            # Assign the first index value
            if myID == 0:
                myIndexList[0] = 0
            else:
                myIndexList[0] = sum(numDOF_all[:myID])

            # Assign the last index value
            myIndexList[1] = myIndexList[0] + numDOF

            # Now offset the indices to take into acount the multiple pyWarp instances
            myIndexList[0] = myIndexList[0] + instanceOffset
            myIndexList[1] = myIndexList[1] + instanceOffset

            # Store this interval in the interval list
            indexList.append(myIndexList)

            # Increment the offset for the next instance
            instanceOffset = instanceOffset + np.sum(numDOF_all)

    #--------------------------------------------

    # Distribute values to the global vector
    for instanceID in range(len(indexList)):

        # Get index interval
        myIndexList = indexList[instanceID]

        # Every proc will send its instance slice to the global vector
        # Citing the PETSc example 2f:
        # "Each processor can contribute any vector entries,
        #  regardless of which processor "owns" them; any nonlocal
        #  contributions will be transferred to the appropriate processor
        #  during the assembly process."
        globalVec.setValues(range(myIndexList[0],myIndexList[1]),
                            localList[instanceID])

        # Assemble vector, using the 2-step process:
        #  VecAssemblyBegin(), VecAssemblyEnd()
        # Computations can be done while messages are in transition
        # by placing code between these two statements.
        # But we will not do this here...
        globalVec.assemblyBegin()
        globalVec.assemblyEnd()

    #--------------------------------------------
    # RETURNS
    return indexList
