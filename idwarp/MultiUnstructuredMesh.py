#!/usr/bin/python
"""

The MultiUnstructuredMesh module is used for creating multiple
USMesh instances, typically from a structured overset mesh.

It contains the following classes:

MultiUSMesh: General class for working with multi-component meshes

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
# =============================================================================
# Imports
# =============================================================================
import os
import numpy as np
from random import randint
from mpi4py import MPI
from .MExt import MExt
from . import USMesh, USMesh_C
from baseclasses.utils import Error

try:
    from cgnsutilities import cgnsutilities as cs
except ImportError:
    cs = None


# =============================================================================
# MultiUnstructuredMesh class
# =============================================================================


class MultiUSMesh(object):
    """
    This mesh object is designed to support independent deformation of
    multiple overset component meshes.
    """

    def __init__(self, CGNSFile, optionsDict, comm=None, dtype="d", debug=False):
        """
        Create the MultiUSMesh object.

        INPUTS:

        CGNSFile: string -> file name of the CGNS file. This CGNS file should be generated with
        cgns_utils combine, so that the domain names have the appropriate convention. That is,
        domains will have the same name as their original files. Domains that share the same name
        will be grouped to make an IDWarp instance.

        optionsDict: dictionary of dictionaries -> Dictionary containing dictionaries that will
        be used to initialize multiple IDWarp instances. The keys are domain names and the
        values are dictionaries of standard IDWarp options that will be applied to this domain.
        The domains of the full CGNS file that do not have a corresponding entry in optionsDict will
        not be warped. For instance, if the CGNS file has the domains wing.00000, wing.00001, and wing.00002
        associated with a wing mesh that we want to warp, then optionsDict should have an entry for 'wing'.

        Ney Secco 2017-02
        """

        # Check if cs was imported correctly:
        if cs is None:
            raise Error("cgns_utils could not be loaded correctly. MultiUSMesh " "requires cgns_utils to function.")

        # Assign communicator if we do not have one yet
        if comm is None:
            comm = MPI.COMM_WORLD

        # Check if warp has already been set by the complex version
        try:
            self.warp
        except AttributeError:
            curDir = os.path.basename(os.path.dirname(os.path.realpath(__file__)))
            self.warp = MExt("libidwarp", curDir, debug=debug)._module

        # Store communicator
        self.comm = comm

        # Store original file name
        self.CGNSFile = CGNSFile

        # Get rank of the current proc
        self.myID = self.comm.Get_rank()

        # Store scalar type
        self.dtype = dtype

        # Set a random prefix to avoid I/O clashes between instances
        prefix = f"tmp{randint(1, 1e5)}"

        # Only the root processor will take the combined CGNS file
        # and explode it by instance.
        if self.myID == 0:

            # Initialize list to store the block IDs that belong to each IDWarp instance.
            # For example, suppose that our combined CGNS file has 21 blocks.
            # Blocks 1 to 5 belong to the fuselage
            # Blocks 6 to 12 belong to the wing
            # Blocks 13 to 21 belong to the background mesh
            # Then cgnsBlockIntervals = [[0,5],[5,12],[12,21]
            self.cgnsBlockIntervals = []

            # Initialize array to store volume nodes CGNS intervals for each instance
            self.cgnsVolNodeIntervals = []

            # Initialize block counter
            blockCounter = 0

            # Initialize node counter
            nodeCounter = 0

            # Load the CGNS file
            combined_file = cs.readGrid(CGNSFile)

            # Explode the CGNS file by zones (this will only work if the user used cgns_utils combine
            # to create the input CGNS file, since the explosion uses the domain names)
            grids, zoneNames = cs.explodeByZoneName(combined_file)

            # Save temporary grid files with the exploded zones
            for grid, zoneName in zip(grids, zoneNames):
                grid.writeToCGNS(f"{prefix}_{zoneName}.cgns")

                # Store the number of blocks in each zone
                self.cgnsBlockIntervals.append([blockCounter, blockCounter + len(grid.blocks)])

                # Count the number of nodes (here is degrees of freedom or nodes*3
                totalNodes = 0
                for blk in grid.blocks:
                    totalNodes += blk.dims[0] * blk.dims[1] * blk.dims[2] * 3

                # Store the number of volume nodes in each zone
                self.cgnsVolNodeIntervals.append([nodeCounter, nodeCounter + totalNodes])

                # Update block counter
                blockCounter = blockCounter + len(grid.blocks)

                # Update node counter
                nodeCounter = nodeCounter + totalNodes

            # Delete grids to free space
            del grids
            del combined_file

        else:

            # Initialize variables to get results in the end
            zoneNames = None
            self.cgnsBlockIntervals = None
            self.cgnsVolNodeIntervals = None

        # Send information to all procs
        zoneNames = self.comm.bcast(zoneNames, root=0)
        self.cgnsBlockIntervals = self.comm.bcast(self.cgnsBlockIntervals, root=0)
        self.cgnsVolNodeIntervals = self.comm.bcast(self.cgnsVolNodeIntervals, root=0)

        # Get names for nearfield meshes.
        # The nearfield mesh names will be the keys of the options dictionary.
        nearfieldNames = optionsDict.keys()

        # Initialize list of IDWarp instances
        self.meshes = []

        # Initialize list to hold indices of the background zones
        self.backgroundInstanceIDs = []

        # Loop over all zones that we found in the combined CGNS file
        for zoneNumber, zoneName in enumerate(zoneNames):

            # Check if the zone belongs to a nearfield mesh
            if zoneName in nearfieldNames:

                # ------------------------------------------------------
                # READING NEARFIELD MESHES (The ones that will be warped)
                #

                # Assign the name of the temporary CGNS file to the options.
                # This is the file that contains the mesh o a single component.
                # Remember that we should use the temporary grid file.
                optionsDict[zoneName]["gridFile"] = f"{prefix}_{zoneName}.cgns"

                # Initialize an IDWarp instance with the current options
                if self.dtype == "d":
                    currMesh = USMesh(options=optionsDict[zoneName], comm=self.comm)
                elif self.dtype == "D":
                    currMesh = USMesh_C(options=optionsDict[zoneName], comm=self.comm)

            else:

                # We have a background mesh

                # Regenerate the temporary filename for the background grid
                bgFile = f"{prefix}_{zoneName}.cgns"

                # ------------------------------------------------------
                # READING BACKGROUND MESHES

                # =========================================================#
                # THIS IS A MESSY (HOPEFULLY TEMPORARY) WAY OF LOADING THE
                # BACKGROUND MESH NODES. IF YOU COME UP WITH A BETTER WAY
                # TO GET volNodes, PLEASE ADD IT HERE.
                # volNodes is a flattened vector that contains the background
                # mesh volume nodes that belong to the current proc.

                # Let's try using IDWarp's CGNS loader to extract the bakground nodes.
                # However, we will have to trick IDWarp in order to do this, since it
                # expects a surface mesh in the file.
                # So we will make a copy of the background mesh file, assign an arbitrary
                # wall surface, and then load it with IDWarp

                # Only the root proc will modify the input file
                if self.myID == 0:

                    # Make a copy of the background mesh file
                    os.system(f"cp {bgFile} {prefix}_bg_file.cgns")

                    # Create a temporary BC file
                    with open(f"{prefix}_bcdata.dat", "w") as fid:
                        fid.write("1 iLow BCwall wall\n")

                    # Use CGNS utils to modify the BCs
                    os.system(f"cgns_utils overwriteBC {prefix}_bg_file.cgns {prefix}_bcdata.dat")

                # Create dummy set of options just to load the CGNS file
                dummyOptions = {
                    "gridFile": f"{prefix}_bg_file.cgns",
                }

                # Initialize an IDWarp instance with the current options
                if self.dtype == "d":
                    currMesh = USMesh(options=dummyOptions, comm=self.comm)
                elif self.dtype == "D":
                    currMesh = USMesh_C(options=dummyOptions, comm=self.comm)

                # Initialize a dummy surface in the background mesh
                """
                if self.myID == 0:
                    print('===========================================')
                    print('ATTENTION: This is a dummy initialization for background mesh warping.')
                pts = np.array([[1.0, 0.0, 1.0],
                                [2.0, 0.0, 1.0],
                                [2.0, 1.0, 1.0],
                                [1.0, 1.0, 1.0]])*(self.myID+1)
                conn = np.array([0,1,2,3])
                faceSizes = np.array([4])
                currMesh.setSurfaceDefinition(pts, conn, faceSizes)
                if self.myID == 0:
                    print('Dummy initialization is Done!')
                    print('===========================================')
                """
                if self.myID == 0:
                    print("===========================================")
                    print("ATTENTION: This is a dummy initialization for background mesh warping.")

                currMesh._setInternalSurface()

                if self.myID == 0:
                    print("Dummy initialization is Done!")
                    print("===========================================")

                # The root proc can remove the temporary files
                if self.myID == 0:

                    # Make a copy of the background mesh file
                    os.system(f"rm {prefix}_bg_file.cgns")
                    os.system(f"rm {prefix}_bcdata.dat")

                # Store the ID of this zone
                self.backgroundInstanceIDs.append(zoneNumber)

            # Append the instance to the list.
            # We will store even the background mesh instances for now,
            # but we will delete them as soon as we call self.setExternalMeshIndices().
            self.meshes.append(currMesh)

        # Now the root proc can remove the temporary grid files
        if self.myID == 0:
            for zoneName in zoneNames:
                os.system(f"rm {prefix}_{zoneName}.cgns")

        # ------------------------------------------------------
        # Initialize other fields for completness
        self.numSurfNodes = None  # How many solver surface nodes we have in the current proc, for all instances
        self.numVolNodes = None  # How many solver volume nodes we have in the current proc, for all instances
        self.cgnsVolNodeMasks = []  # Mask used to filter which volume nodes given by the solver belong to each instance

    def getSurfaceCoordinates(self):
        """Returns all defined surface coordinates on this processor, with the Solver ordering

        Returns
        -------
        pts : numpy array size (N,3)
            Specified surface coordinates residing on this
            processor. This may be empty array, size (0,3)

        Ney Secco 2017-02
        """

        # Check if the user already set the surface definition
        if self.numSurfNodes is None:
            if self.comm.Get_rank() == 0:
                raise NameError("IDWarpMulti surface is not initialized. Run self.setSurfaceDefinition first.")

        # Initialize the array of points
        pts = np.zeros((self.numSurfNodes, 3), dtype=self.dtype)

        # Loop over every instance to get their contributions
        for instanceID, mesh in enumerate(self.meshes):

            # Get current set of points
            currPts = mesh.getSurfaceCoordinates()

            # Assign this set of points to the full one
            pts[self.filtered2fullMaps[instanceID], :] = currPts

        # Return set of points
        return pts

    def setSurfaceCoordinates(self, pts):
        """Sets all surface coordinates on this processor, with pts given in solver ordering

        Parameters
        ----------
        pts : numpy array, size(N, 3)
            The coordinate to set. This MUST be exactly the same size as
            the array obtained from getSurfaceCoordinates()

        Ney Secco 2017-02
        """

        # Check if the user already set the surface definition
        if self.numSurfNodes is None:
            if self.comm.Get_rank() == 0:
                raise NameError("IDWarpMulti surface is not initialized. Run self.setSurfaceDefinition first.")

        # Loop over every mesh object to get a slice of the surface point array
        for instanceID, mesh in enumerate(self.meshes):

            # Extract set of points that belong to this instance
            currPts = pts[self.filtered2fullMaps[instanceID], :]

            # Set points to the current instance
            mesh.setSurfaceCoordinates(currPts)

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

        # Here is the main logic of the function:
        # - We have all volume node IDs that belong to each instance, in CGNS ordering,
        #   stored in self.cgnsVolNodeIntervals. All procs have this same array.
        # - ADflow gives us the CGNS indices that the current proc needs (this is ind).
        #   Thus we can use self.cgnsVolNodeIntervals to determine which instance is responsible by
        #   each ind entry.
        # - So we create masks so that we can take pieces from the full volume vector and
        #   send them to each instance.

        # Print log
        if self.myID == 0:
            print("")
            print("Mapping solver volume nodes to IDWarpMulti volume nodes")

        # Save the number of volume nodes
        self.numVolNodes = len(ind)

        # Here we will initialize a vector of volume nodes.
        # We will populate this volume node vector with the current mesh state, including
        # the background mesh nodes. We will use this vector to initialize any new volume
        # node vector returned by IDWarpMulti when we call self.getSolverGrid().
        # This will allow us to delete all IDWarp instances related to the background nodes
        # since they will remain unchanged.
        self.defaultSolverGrid = np.zeros(self.numVolNodes, dtype=self.dtype)

        # Loop over every instance
        for instanceID, mesh in enumerate(self.meshes):

            # Get CGNS bounds of the current instance
            startIndex = self.cgnsVolNodeIntervals[instanceID][0]
            endIndex = self.cgnsVolNodeIntervals[instanceID][1]

            # Create and store mask that will select the elements of the volume array that belong
            # to the current instance
            currMask = (ind >= startIndex) * (ind < endIndex)

            # Store the mask
            self.cgnsVolNodeMasks.append(currMask)

            # Gather elements of ind are within the CGNS interval owned by this instance
            ownedInd = ind[currMask]

            # Offset the indices, since the instance CGNS ordering starts at 0 since it does not
            # know about the other instances
            ownedInd = ownedInd - startIndex

            # Now we can set the external mesh indices in the current instance
            mesh.setExternalMeshIndices(ownedInd)

            # Retrieve background mesh nodes to start populating the default volume node vector
            if instanceID in self.backgroundInstanceIDs:
                self.defaultSolverGrid[currMask] = mesh.getSolverGrid()

        # Now that we populated the default volume node vector, we can delete the background
        # mesh instances
        for instanceID in reversed(range(len(self.meshes))):

            # Check if the current instance is a background mesh
            if instanceID in self.backgroundInstanceIDs:

                # Delete current instance
                del self.meshes[instanceID]

                # Delete it CGNS intervals
                del self.cgnsVolNodeMasks[instanceID]
                del self.cgnsVolNodeIntervals[instanceID]
                del self.cgnsBlockIntervals[instanceID]

        # State that we removed all background meshes
        self.backgroundInstanceIDs = None

        # Print log
        if self.myID == 0:
            print(" Done!")
            print("")

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

        Ney Secco 2017-02
        """

        # Initialize array that will hold all volume nodes in the current proc in solver ordering.
        # We use a copy of the default solver grid vector, which contains the initial mesh (including the BG meshes)
        # as a baseline. Later we will update just the nearfield nodes.
        solverGrid = self.defaultSolverGrid.copy()

        # Loop over each instance to gather volume nodes
        for instanceID, mesh in enumerate(self.meshes):

            # Get volume nodes of the current instance
            currVolNodes = mesh.getSolverGrid()

            # Assign volume nodes to the full array using the mask of the current instance
            solverGrid[self.cgnsVolNodeMasks[instanceID]] = currVolNodes

        return solverGrid

    def setSurfaceDefinition(self, pts, conn=None, faceSizes=None, cgnsBlockIDs=None, distTol=1e-6):
        """This is the master function that determines the definition of the
        surface to be used for the mesh movement. This surface may be
        supplied from an external solver (such as ADflow) or it may be
        generated by IDWarp internally.

        Parameters
        ----------
        pts : array, size (M, 3)
            Nodes on this processor
        conn : int array, size  (sum(faceSizes))
            Connectivity of the nodes on this processor
        faceSizes : int Array size (N)
            Treat the conn array as a flat list with the faceSizes giving
            connectivity offset for each element.
        cgnsBlockIDs : int Array size (N)
            Block ID, in CGNS ordering, that contains each element.
        distTol: Distance tolerance to flag that a given surface node does not
                 belong to the current IDWarp surface definition in the current proc.

        Ney Secco 2017-02
        """

        # Here is the main logic of this function:
        # - Each proc receives a set of surface patch from ADflow
        # - cgnsBlockIDs indicates which CGNS block contains each surface patch
        # - We stored the CGNS block IDs that belong to each IDWarp instance at the
        #   initialization method.
        # - So we can use the CGNS block ID to separate ADflow surface patches per
        #   IDWarp instance.

        # Print log
        if self.myID == 0:
            print("")
            print("Mapping solver surface nodes to IDWarp volume nodes")

        # Check if the user executed self.setExternalMeshIndices first to remove background meshes
        if self.backgroundInstanceIDs is not None:
            raise NameError("The user should run self.setExternalMeshIndices before self.setSurfaceDefinition.")

        ### Initial steps

        # Get the number of elements
        numElems = len(faceSizes)

        # Check if we have the same number of CGNS indices
        if len(cgnsBlockIDs) != numElems:
            raise ValueError("The user should provide block IDs for every element.")

        # Save the total number of surface nodes
        self.numSurfNodes = pts.shape[0]

        # Make sure we have flat connectivities
        conn = conn.flatten()

        ### We will take the cgnsBlockIDs given by ADflow and convert them into instance IDs.

        # Initialize array to store ID of the instance that has each cell
        self.instanceIDs = np.ones(numElems) * (-1)

        # Now we use the block ID bounds of each instance to flag the instance IDs
        for ii in range(len(self.meshes)):

            # Get indices of the first and last CGNS block that belongs to the current IDWarp instance
            indexStart = self.cgnsBlockIntervals[ii][0]
            indexEnd = self.cgnsBlockIntervals[ii][1]

            # Check which indices are defined within the range of cgns blocks.
            # Then we assign the current instance ID to the correponding cells
            self.instanceIDs[(cgnsBlockIDs >= indexStart) * (cgnsBlockIDs < indexEnd)] = ii

        ### We will send the connectivity subsets to each instance

        # We will need to create subgroups of points and connectivities to each instance since
        # IDWarp cannot receive points that are not connected to any element.

        # Initialize lists that will map the full set of nodes to each subgroup of nodes
        self.filtered2fullMaps = []

        for ii in range(len(self.meshes)):

            # Initialize list to hold connectivities and face sizes for the current instance
            currConn = []
            currFaceSizes = []

            # Initialize counter to track our position along the conn array
            connPos = 0

            # Now loop over every element to gather just the ones that belong to the current proc
            for elemID in range(numElems):

                # Get the size of the current element
                elemFaceSize = faceSizes[elemID]

                # Check if the current element belongs to the current IDWarp instance
                if self.instanceIDs[elemID] == ii:

                    # Store its faceSize
                    currFaceSizes.append(elemFaceSize)

                    # Store its connectivity
                    currConn = currConn + conn[connPos : connPos + elemFaceSize].tolist()

                # Increment position counter for next iteration
                connPos = connPos + elemFaceSize

            # Transform connectivities back to numpy array
            currConn = np.array(currConn, dtype="intc").flatten()
            currFaceSizes = np.array(currFaceSizes, dtype="intc")

            # Now we will use the filtered connectivities to gather only the points
            # that belong to the current instance

            # Get a list of points used by the filtered connectivities. We can remove
            # all repetitions in this case (this is what set() does).
            # Note that the result is equivalent to the mapping between the filtered nodes
            # and the full set of nodes.
            # For instance if filtered2fullMap[ii] = jj, this means that pts[jj,:] = currPts[ii,:]
            filtered2fullMap = list(set(currConn.tolist()))

            # Now we can extract the filtered points from the full array
            currPts = pts[filtered2fullMap, :]

            # Currently, the indices in currConn refer to the full set of points. We need to
            # update these connectivities so that the point to the filtered set of points (currPts).

            # Loop over every mapping:
            for filteredID, fullID in enumerate(filtered2fullMap):

                # Find all elements of currConn that share the same fullID and replace them
                # with the new ID
                currConn[currConn == fullID] = filteredID

            # Now we finally can set surface definition of the current IDWarp instance
            self.meshes[ii].setSurfaceDefinition(currPts, currConn, currFaceSizes)

            # Save the mapping for future uses
            self.filtered2fullMaps.append(filtered2fullMap)

        # Print log
        if self.myID == 0:
            print("Done")
            print("")

    def getWarpGrid(self):
        """
        Return the current grids. This function is typically unused. See
        getSolverGrid for the more useful interface functionality.

        This only returns the nearfield meshes.

        Returns
        -------
        volNodesList, list of 1D numpy arrays, real: These are the local volume nodes (in a flat 1D array)
        of each instance. That is, volNodesList[i] has the volume nodes stored in the local proc for
        the i-th IDWarp instance.

        numCoorTotal: The total number of coordinates, across all procs and IDWarp instances.

        Ney Secco 2017-02
        """

        # Initialize list to hold volume nodes (in the current proc) of all instances.
        # That is, volNodesList[i] gives the volume nodes of the i-th IDWarp instance that belong
        # to the current proc.
        volNodesList = []

        # Initialize counter to store the total number of coordinates (numCoor = 3*numNodes)
        # of all volume mesh.
        numCoorTotal = 0

        # Loop over the multiple CGNS files to initialize the corresponding IDWarp instances
        for currMesh in self.meshes:

            # Get volume nodes.
            # volNodes is a flattened vector that contains the background
            # mesh volume nodes that belong to the current proc.
            volNodes = currMesh.getCommonGrid()

            # Store the nodes of the current instance in the list
            volNodesList.append(volNodes)

            # Get number of coordinates on the current processor, for the current IDWarp instance.
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
        mesh-warp derivative computation. Note that the same steps
        used in this function are done at warpDeriv, so in theory you
        could save time by just saving the output of warpDeriv.

        Here we accumulate all seeds coming from each IDWarp instance

        Returns
        -------
        dXs :  numpy array
            The specific components of dXs. size(N,3). This the same
            size as the array obtained with getSurfaceCoordiantes(). N may
            be zero if this processor does not have any surface coordinates.
            These can be, for instance, reverse AD seeds.

        Ney Secco 2017-03
        """

        # Initialize array to hold seeds of all surface point of this proc
        dXs = np.zeros((self.numSurfNodes, 3), dtype=self.dtype)

        # Loop over all instances
        for instanceID, mesh in enumerate(self.meshes):

            # Get the surface seeds of the current instance
            curr_dXs = mesh.getdXs()

            # Assign seeds of the current instance to the general array
            dXs[self.filtered2fullMaps[instanceID]] = dXs[self.filtered2fullMaps[instanceID]] + curr_dXs

        return dXs

    def warpMesh(self):
        """
        This calls the mesh warping method for each IDWarp instance.

        This will update the volume coordinates internally in each instance.

        Ney Secco 2017-02
        """

        # Print log
        if self.myID == 0:
            print("")
            print("Starting IDWarpMulti mesh warping")

        # Set mesh counter
        meshCounter = 1

        # Loop over all instances
        for mesh in self.meshes:

            # Print log
            if self.myID == 0:
                print("")
                print(" warping mesh", meshCounter, "of", len(self.meshes))

            # Warp current instance
            mesh.warpMesh()

            # Print log
            if self.myID == 0:
                print("")
                print(" Done")

            # Increment counter
            meshCounter = meshCounter + 1

        # Print log
        if self.myID == 0:
            print("")
            print("IDWarpMulti successfully warped all instances!")
            print("")

    def warpDeriv(self, dXv, solverVec=True):
        """Compute the warping derivative (dXv/dXs^T)*Vec (where vec is the
        dXv argument to this function.

        This is the main routine to compute the mesh warping
        derivative.

        Parameters
        ----------
        dXv :  numpy array
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
        dXs :  numpy array
            The specific components of dXs. size(N,3). This the same
            size as the array obtained with getSurfaceCoordiantes(). N may
            be zero if this processor does not have any surface coordinates.
            These can be, for instance, reverse AD seeds.

        Ney Secco 2017-03
        """

        # Print log
        if self.myID == 0:
            print("")
            print("Starting IDWarpMulti reverse AD")

        # ---------------------------------------------------
        # RUN REVERSE AD ON EACH INSTANCE

        # Loop over all instances
        for instanceID, mesh in enumerate(self.meshes):

            # Print log
            if self.myID == 0:
                print("")
                print(" Working on mesh", instanceID + 1, "of", len(self.meshes))

            # Get current seeds
            curr_dXv = dXv[self.cgnsVolNodeMasks[instanceID]]

            # Run reverse AD.
            # This will update the surface seeds inside the mesh object
            mesh.warpDeriv(curr_dXv, solverVec=solverVec)

            # Print log
            if self.myID == 0:
                print("")
                print(" Done")

        # Print log
        if self.myID == 0:
            print("")
            print("IDWarpMulti successfully finished reverse AD on all instances!")
            print("")

        # Get derivative seeds
        dXs = self.getdXs()

        # Return the computed surface seeds
        return dXs

    def warpDerivFwd(self, dXs):
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

        Ney Secco 2017-03
        """

        # Print log
        if self.myID == 0:
            print("")
            print("Starting IDWarpMulti forward AD")

        # ---------------------------------------------------
        # SPLIT SURFACE SEEDS ACROSS ALL INSTANCES AND GATHER
        # ALL VOLUME SEEDS

        # Initialize list to gather volume seeds of all instances
        dXv = np.zeros(self.numVolNodes, dtype=self.dtype)

        # Loop over every mesh object to get a slice of the surface seeds array
        for instanceID, mesh in enumerate(self.meshes):

            # Print log
            if self.myID == 0:
                print("")
                print(" Working on mesh", instanceID + 1, "of", len(self.meshes))

            # Get current surface node forward AD seeds
            curr_dXs = dXs[self.filtered2fullMaps[instanceID]]

            # Run the forward AD code for this instance to get volume node seeds.
            curr_dXvWarp = mesh.warpDerivFwd(curr_dXs)

            # Add seeds to the full volume vector
            dXv[self.cgnsVolNodeMasks[instanceID]] = curr_dXvWarp

            # Print log
            if self.myID == 0:
                print("")
                print(" Done")

        # Print log
        if self.myID == 0:
            print("")
            print("IDWarpMulti successfully finished forward AD on all instances!")
            print("")

        # ---------------------------------------------------
        # RETURNS
        return dXv

    # ==========================================================================
    #                        Output Functionality
    # ==========================================================================
    def writeGrid(self, baseName=None):
        """
        Write the grids of each instance

        Parameters
        ----------
        baseName : str or None
            a base namve that will be used to generate filenames for
            all instance mesh files.

        Ney Secco 2017-03
        """

        # We need to reexplode the original file
        if self.myID == 0:

            # Load the CGNS file
            combined_file = cs.readGrid(self.CGNSFile)

            # Explode the CGNS file by zones (this will only work if the user used cgns_utils combine
            # to create the input CGNS file, since the explosion uses the domain names)
            grids, zoneNames = cs.explodeByZoneName(combined_file)

            # Save temporary grid files with the exploded zones
            for grid, zoneName in zip(grids, zoneNames):
                grid.writeToCGNS("_" + zoneName + ".cgns")

            # Delete grids to free space
            del grids
            del combined_file

        # Assign baseName if user provided none
        if baseName is None:
            baseName = "IDWarpMulti"

        # Loop over all instances to write their meshes
        for instanceID, mesh in enumerate(self.meshes):

            # Generate a fileName
            fileName = baseName + "_inst%03d.cgns" % (instanceID)

            # Call function to write the mesh of the current instance
            mesh.writeGrid(fileName)

        # Now the root proc can remove the temporary grid files
        if self.myID == 0:
            for zoneName in zoneNames:
                os.system("rm _" + zoneName + ".cgns")

    # =========================================================================
    #                     Utility functions
    # =========================================================================

    def __del__(self):
        """Release all the mesh warping memory. This should be called
        automatically when the object is garbage collected."""
        for mesh in self.meshes:
            mesh.warp.releasememory()
