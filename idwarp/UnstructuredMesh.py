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
# =============================================================================
# Imports
# =============================================================================
import os
import shutil
import numpy as np
from mpi4py import MPI
from .MExt import MExt
from baseclasses import BaseSolver
from baseclasses.utils import Error


class Warning(object):
    """
    Format a warning message
    """

    def __init__(self, message):
        msg = "\n+" + "-" * 78 + "+" + "\n" + "| IDWarp Warning: "
        i = 17
        for word in message.split():
            if len(word) + i + 1 > 78:  # Finish line and start new one
                msg += " " * (78 - i) + "|\n| " + word + " "
                i = 1 + len(word) + 1
            else:
                msg += word + " "
                i += len(word) + 1
        msg += " " * (78 - i) + "|\n" + "+" + "-" * 78 + "+" + "\n"
        print(msg)


# =============================================================================
# UnstructuredMesh class
# =============================================================================


class USMesh(BaseSolver):
    """
    This is the main Unstructured Mesh. This mesh object is designed
    to interact with an structured or unstructured CFD solver though a
    variety of interface functions.
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

        name = "IDWarp"
        category = "Volume mesh warping"

        if comm is None:
            comm = MPI.COMM_WORLD

        # Default options for mesh warping
        defOpts = self._getDefaultOptions()

        if options is None:
            raise Error(
                "The 'options' keyword argument is *NOT* "
                "optional. An options dictionary must be passed upon "
                "creation of this object"
            )

        # Initialize the inherited BaseSolver
        super().__init__(name, category, defaultOptions=defOpts, options=options, comm=comm)

        # aExp and bExp are not fully implemented
        if self.getOption("aExp") != 3.0 or self.getOption("bExp") != 5.0:
            raise Error(
                "The aExp and bExp options are currently not implemented "
                "and should not be modified from their default values."
            )

        self.printOptions()

        # Check if warp has already been set if this has been
        # inherited to complex version
        try:
            self.warp
        except AttributeError:
            curDir = os.path.basename(os.path.dirname(os.path.realpath(__file__)))
            self.warp = MExt("libidwarp", curDir, debug=debug)._module

        # Initialize PETSc if not done so
        self.warp.initpetsc(self.comm.py2f())

        # Set realtype of 'd'. 'D' is used in Complex and set in
        # UnstructuredMesh_C.py
        self.dtype = "d"

        # Set Fortran options values
        self._setMeshOptions()

        # Initialize various bits of stored information
        self.OFData = {}
        self.warpInitialized = False
        self.faceSizes = None
        self.fileType = self.getOption("fileType")
        fileName = self.getOption("gridFile")

        # Determine how to read
        if self.fileType == "CGNS":
            # Determine type of CGNS mesh we have
            self.warp.readcgns(fileName)
        elif self.fileType == "OpenFOAM":
            self._readOFGrid(fileName)
        elif self.fileType == "PLOT3D":
            self.warp.readplot3d(fileName)

    @staticmethod
    def _getDefaultOptions():
        defOpts = {
            "gridFile": [(str, type(None)), None],
            "fileType": [str, ["CGNS", "OpenFOAM", "PLOT3D"]],
            "specifiedSurfaces": [(list, type(None)), None],
            "symmetrySurfaces": [(list, type(None)), None],
            "symmetryPlanes": [(list, type(None)), None],
            "aExp": [float, 3.0],
            "bExp": [float, 5.0],
            "LdefFact": [float, 1.0],
            "alpha": [float, 0.25],
            "errTol": [float, 0.0005],
            "evalMode": [str, ["fast", "exact"]],
            "symmTol": [float, 1e-6],
            "useRotations": [bool, True],
            "zeroCornerRotations": [bool, True],
            "cornerAngle": [float, 30.0],
            "restartFile": [(str, type(None)), None],
            "bucketSize": [int, 8],
        }
        return defOpts

    def getSurfaceCoordinates(self):
        """Returns all defined surface coordinates on this processor

        Returns
        -------
        coords : numpy array size (N,3)
            Specified surface coordinates residing on this
            processor. This may be empty array, size (0,3)
        """
        self._setInternalSurface()
        coords = np.zeros((self.nSurf, 3), self.dtype)
        self.warp.getsurfacecoordinates(np.ravel(coords))

        return coords

    def setSurfaceCoordinates(self, coordinates):
        """Sets all surface coordinates on this processor

        Parameters
        ----------
        coordinates : numpy array, size(N, 3)
            The coordinate to set. This MUST be exactly the same size as
            the array obtained from getSurfaceCoordinates()
        """
        if len(coordinates) != self.nSurf:
            raise Error(
                "Incorrect length of coordinates supplied to "
                "setSurfaceCoordinates on proc %d. Expected "
                "aray of length %d, received an array of length "
                "%d." % (self.comm.rank, self.nSurf, len(coordinates))
            )

        self.warp.setsurfacecoordinates(np.ravel(coordinates))

    def setExternalMeshIndices(self, ind):
        """
        Set the indicies defining the transformation of an external
        solver grid to the original CGNS grid. This is required to use
        USMesh functions that involve the word "Solver" and warpDeriv.
        The indices must be zero-based.

        Parameters
        ----------
        ind : numpy integer array
            The list of indicies this processor needs from the common mesh file
        """

        self.warp.setexternalmeshindices(ind)

    def setSurfaceDefinition(self, pts, conn, faceSizes, cgnsBlockID=None):
        """This is the master function that determines the definition of the
        surface to be used for the mesh movement. This surface may be
        supplied from an external solver (such as SUmb) or it may be
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
        cgnsBlockID : dummy argument.
            This argument is not used at all. It is here just to have the
            same API as the IDWarpMulti class.
        """

        conn = np.array(conn).flatten()

        # Call the fortran initialze warping command with  the
        # coordinates and the patch connectivity given to us.
        if self.getOption("restartFile") is None:
            restartFile = ""
        else:
            restartFile = self.getOption("restartFile")

        # The symmetry conditions but be set before initializing the
        # warping.
        self._setSymmetryConditions()

        allPts = np.vstack(self.comm.allgather(pts))
        allFaceSizes = np.hstack(self.comm.allgather(faceSizes))
        # Be careful with conn...need to increment by the offset from
        # the points.
        ptSizes = self.comm.allgather(len(pts))
        offsets = np.zeros(len(ptSizes), "intc")
        offsets[1:] = np.cumsum(ptSizes)[:-1]
        allConn = np.hstack(self.comm.allgather(conn + offsets[self.comm.rank]))

        # Next point reduce all nodes:
        uniquePts, link, nUnique = self.warp.pointreduce(allPts.T, 1e-12)
        uniquePts = uniquePts.T[0:nUnique]

        # Now sort them by z, y and x to ensure that the tree is
        # always created the same way. It is currently sensitivte to
        # the ordering of the original nodes. Ie, if you reorder the
        # nodes you get a differnet tree which gives slightly
        # differnet warp values/derivatives.

        ind = np.lexsort(uniquePts.T)
        uniquePts = uniquePts[ind]

        # Update the link array such that it points back to the original list.
        inv = np.zeros_like(ind, "intc")
        rng = np.arange(len(ind)).astype("intc")

        inv[ind[rng]] = rng
        link = inv[link - 1] + 1
        self.warp.initializewarping(np.ravel(pts), uniquePts.T, link, allFaceSizes, allConn, restartFile)
        self.nSurf = len(pts)
        self.warpInitialized = True

        # Clean-up patch data no longer necessary.
        self.warp.deallocatepatchio()

    def setSurfaceDefinitionFromFile(self, surfFile):
        """Set the surface definition for the warping from a
        multiblock PLOT3D surface file

        Parameters
        ----------
        surfFile : filename of multiblock PLOT3D surface file.
        """

        # Read the PLOT3D surface file
        self.warp.readplot3dsurface(surfFile)

        if self.comm.rank == 0:
            # Extract out the data we want from the module. Note that
            # becuase we are in python, convert conn to 0 based.
            pts = self.warp.plot3dsurface.pts.T.copy()
            conn = self.warp.plot3dsurface.conn.T.copy() - 1
            faceSizes = 4 * np.ones(len(conn))
        else:
            pts = np.zeros((0, 3))
            conn = np.zeros((0, 4))
            faceSizes = 4 * np.ones(0)

        # Run the regular setSurfaceDefinition routine. Note that only
        # the root proc has non-zero length arrays. That's ok.
        self.setSurfaceDefinition(pts, conn, faceSizes)

        # Cleanup the fortran memory we used to read the surface
        self.warp.plot3dsurface.pts = None
        self.warp.plot3dsurface.conn = None

    def setSurfaceFromFile(self, surfFile):
        """Update the internal surface surface coordinates using an
        external PLOT3D surface file. This can be used in an analogous
        way to setSurfaceDefinitionFromFile. The 'sense' of the file
        must be same as the file used with
        setSurfaceDefinitionFromFile. That means, the same number of
        blocks, with the same sizes, in the same order.

        Parameters
        ----------
        surfFile: filename of multiblock PLOT3D surface file'
        """

        # Read the PLOT3D surface file
        self.warp.readplot3dsurface(surfFile)

        if self.comm.rank == 0:
            # Extract out the data we want from the module
            pts = self.warp.plot3dsurface.pts.T.copy()
        else:
            pts = np.zeros((0, 3))

        # Do the regular update
        self.setSurfaceCoordinates(pts)

        # Cleanup the memory we used to read the surface
        self.warp.plot3dsurface.pts = None
        self.warp.plot3dsurface.conn = None

    def getSolverGrid(self):
        """Return the current grid in the order specified by
        setExternalMeshIndices(). This is the main routine for
        returning the defomed mesh to the external CFD solver.

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
        """Return the grid in the original ordering. This is required for the
        OpenFOAM tecplot writer since the connectivity is only known
        in this ordering.
        """
        return self.warp.getcommonvolumecoordinates(self.warp.griddata.commonmeshdof)

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
        dXs = np.zeros((self.nSurf, 3), self.dtype)
        self.warp.getdxs(np.ravel(dXs))

        return dXs

    def warpMesh(self):
        """
        Run the applicable mesh warping strategy.

        This will update the volume coordinates to match surface
        coordinates set with setSurfaceCoordinates()
        """
        # Initialize internal surface if one isn't already set.
        self._setInternalSurface()

        # Set Options
        self._setMeshOptions()

        # Now run mesh warp command
        self.warp.warpmesh()

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

        dXvWarp = np.zeros(self.warp.griddata.warpmeshdof, dtype=self.dtype)
        self.warpMesh()
        self.warp.warpderivfwd(np.ravel(dXs), dXvWarp)
        if solverVec:
            dXv = np.zeros(self.warp.griddata.solvermeshdof, dtype=self.dtype)
            self.warp.warp_to_solver_grid(dXvWarp, dXv)
            return dXv
        else:
            return dXvWarp

    def verifyWarpDeriv(self, dXv=None, solverVec=True, dofStart=0, dofEnd=10, h=1e-6, randomSeed=314):
        """Run an internal verification of the solid warping
        derivatives"""

        if dXv is None:
            np.random.seed(randomSeed)  # 'Random' seed to ensure runs are same
            dXvWarp = np.random.random(self.warp.griddata.warpmeshdof, dtype=self.dtype)
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

            Filename for grid. Should end in .cgns for CGNS files. For
            PLOT3D whatever you want. It is not optional for
            CGNS/PLOT3D. It is not required for OpenFOAM meshes. This
            call will update the 'points' file.
        """
        if self.fileType in ["CGNS", "PLOT3D"]:
            if fileName is None:
                raise Error("fileName is not optional for writeGrid with " "gridType of CGNS or PLOT3D")

            if self.fileType == "CGNS":
                # Copy the default and then write
                if self.comm.rank == 0:
                    shutil.copy(self.getOption("gridFile"), fileName)
                self.warp.writecgns(fileName)

            elif self.fileType == "PLOT3D":
                self.warp.writeplot3d(fileName)

        elif self.fileType == "OpenFOAM":
            self._writeOpenFOAMVolumePoints(self.getCommonGrid())

    def writeOFGridTecplot(self, fileName):
        """
        Write the current OpenFOAM grid to a Tecplot FE polyhedron file.
        This is generally used for debugging/visualization purposes.

        Parameters
        ----------
        fileName : str
            Filename to use. Should end in .dat for tecplot ascii file.
        """
        if not self.OFData:
            Warning("Cannot write OpenFOAM tecplot file since there is " "no OpenFOAM data present")
            return

        if self.comm.size == 1:
            f = open(fileName, "w")
        else:
            fname, fext = os.path.splitext(fileName)
            fileName = fname + "%d" % self.comm.rank + fext
            f = open(fileName, "w")

        # Extract the data we need from the OFDict to make the code a
        # little easier to read:
        nodes = self.getCommonGrid()
        nodes = nodes.reshape((len(nodes) // 3, 3))
        nPoints = len(nodes)

        faces = self.OFData["faces"]
        nFaces = len(faces)

        owner = self.OFData["owner"]
        nCells = np.max(owner) + 1

        neighbour = self.OFData["neighbour"]
        nNeighbours = len(neighbour)

        # write the node numbers for each face
        faceNodeSum = 0
        for i in range(nFaces):
            faceNodeSum += len(faces[i])

        # Tecplot Header Info:
        f.write('TITLE = "Example Grid File"\n')
        f.write("FILETYPE = GRID\n")
        f.write('VARIABLES = "X" "Y" "Z"\n')
        f.write("Zone\n")
        f.write("ZoneType=FEPOLYHEDRON\n")
        f.write("NODES=%d\n" % nPoints)
        f.write("FACES=%d\n" % nFaces)
        f.write("ELEMENTS=%d\n" % nCells)
        f.write("TotalNumFaceNodes=%d\n" % faceNodeSum)
        f.write("NumConnectedBoundaryFaces=%d\n" % 0)
        f.write("TotalNumBoundaryConnections=%d\n" % 0)

        # Write the points to file in block data format
        for idim in range(3):
            for j in range(nPoints):
                f.write("%g\n" % nodes[j, idim])

        # Write the number of nodes in each face
        for i in range(nFaces):
            f.write("%d " % len(faces[i]))
            if i % 300 == 0:
                f.write("\n")

        f.write("\n")
        # write the node numbers for each face (+1 for tecplot 1-based
        # ordering)
        for i in range(nFaces):
            nPointsFace = len(faces[i])
            for j in range(nPointsFace):
                f.write("%d " % (faces[i][j] + 1))
            f.write("\n")

        # write left elements (+1 for 1 based ordering)
        for i in range(nFaces):
            f.write("%d\n" % (owner[i] + 1))

        # write right elements (+1 for 1 based ordering)
        for i in range(nFaces):
            if i < nNeighbours:
                f.write("%d\n" % (neighbour[i] + 1))
            else:
                f.write("%d\n" % (0))

        f.close()

    # =========================================================================
    #                     Internal Private Functions
    # =========================================================================
    def _setInternalSurface(self):
        """This function is used by default if setSurfaceDefinition() is not
        set BEFORE an operation is requested that requires this
        information.
        """

        if self.warpInitialized:
            return

        if self.comm.rank == 0:
            Warning(
                "Using internally generated IDWarp surfaces. If "
                "this mesh object is to be used with an "
                "external solver, ensure the mesh object is "
                "passed to the solver immediatedly after it is created. "
                "The external solver must then call "
                "'setExternalMeshIndices()' and 'setSurfaceDefinition()' "
                "routines."
            )

        conn = []
        pts = []
        faceSizes = []

        if self.fileType == "CGNS":
            if self.comm.rank == 0:

                # Do the necessary fortran preprocessing
                if self.warp.cgnsgrid.cgnsstructured:
                    self.warp.processstructuredpatches()
                else:
                    self.warp.processunstructuredpatches()

                fullConn = self.warp.cgnsgrid.surfaceconn - 1
                fullPts = self.warp.cgnsgrid.surfacepoints

                nPatch = self.warp.cgnsgrid.getnpatch()
                fullPatchNames = []
                for i in range(nPatch):
                    fullPatchNames.append(self.warp.cgnsgrid.getsurf(i + 1).strip().lower())
                # We now have all surfaces belonging to
                # boundary conditions. We need to decide which
                # ones to use depending on what the user has
                # told us.
                surfaceFamilies = set()
                if self.getOption("specifiedSurfaces") is None:
                    # Use all wall surfaces:
                    for i in range(len(self.warp.cgnsgrid.surfaceiswall)):
                        if self.warp.cgnsgrid.surfaceiswall[i]:
                            surfaceFamilies.add(fullPatchNames[i].lower())
                else:
                    # The user has supplied a list of surface families
                    for name in self.getOption("specifiedSurfaces"):
                        surfaceFamilies.add(name.lower())

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
                        curNodeSize = fullPatchSizes[i][0] * fullPatchSizes[i][1]
                        curCellSize = (fullPatchSizes[i][0] - 1) * (fullPatchSizes[i][1] - 1)

                        if fullPatchNames[i] in surfaceFamilies:
                            # Keep track of the families we've actually used
                            usedFams.add(fullPatchNames[i])
                            pts.extend(fullPts[curNodeIndex : curNodeIndex + curNodeSize * 3])
                            conn.extend(fullConn[curCellIndex : curCellIndex + curCellSize * 4] - curOffset)
                        else:
                            # If we skipped, we increment the offset
                            curOffset += curNodeSize

                        curNodeIndex += curNodeSize * 3
                        curCellIndex += curCellSize * 4

                # end for (root proc)

                # Run the common surface definition routine
                pts = np.array(pts).reshape((len(pts) // 3, 3))
                faceSizes = 4 * np.ones(len(conn) // 4, "intc")
                self.setSurfaceDefinition(pts=pts, conn=np.array(conn, "intc"), faceSizes=faceSizes)

            else:  # unstructured
                if self.comm.rank == 0:

                    # Pull out data and convert to 0-based ordering
                    fullPtr = self.warp.cgnsgrid.surfaceptr - 1
                    fullPatchPtr = self.warp.cgnsgrid.surfacepatchptr - 1
                    fullFaceSizes = fullPtr[1:-1] - fullPtr[0:-2]

                    # Now we loop back over the "full" versions of
                    # thing and just take the ones that correspond
                    # to the families we are using.
                    curOffset = 0
                    for i in range(len(fullPatchNames)):

                        # Start/end indices into fullPtr array
                        iStart = fullPatchPtr[i]
                        iEnd = fullPatchPtr[i + 1]

                        # Start/end indices into the fullConn/fullPts array
                        iStart2 = fullPtr[iStart]
                        iEnd2 = fullPtr[iEnd]

                        if fullPatchNames[i] in surfaceFamilies:
                            # Keep track of the families we've actually used
                            usedFams.add(fullPatchNames[i])
                            faceSizes.extend(fullFaceSizes[iStart:iEnd])
                            conn.extend(fullConn[iStart2:iEnd2] - curOffset)
                            pts.extend(fullPts[iStart2 * 3 : iEnd2 * 3])
                        else:
                            # If we skipped, we increment the offset
                            curOffset += iEnd2 - iStart2

                pts = np.array(pts).reshape((len(pts) // 3, 3))
                # Run the common surface definition routine
                self.setSurfaceDefinition(pts=pts, conn=np.array(conn, "intc"), faceSizes=faceSizes)

            # Check if all supplied family names were actually
            # used. The user probably wants to know if a family
            # name was specified incorrectly.
            if self.comm.rank == 0:
                if usedFams < surfaceFamilies:
                    missing = list(surfaceFamilies - usedFams)
                    Warning(
                        "Not all specified surface families that "
                        "were given were found the CGNS file. "
                        "The families not found are %s." % (repr(missing))
                    )
                if len(usedFams) == 0:
                    raise Error(
                        "No surfaces were selected. Check the names "
                        "given in the 'specifiedSurface' option. The "
                        "list of families is %s." % (repr(list(fullPatchNames)))
                    )

        elif self.fileType == "OpenFOAM":

            faceSizes, conn, pts = self._computeOFConn()

            # Run the "external" command
            self.setSurfaceDefinition(pts=pts, conn=conn, faceSizes=faceSizes)

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

        symmList = self.getOption("symmetryPlanes")
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
            if self.fileType == "CGNS":
                if self.comm.rank == 0:

                    # Do the necessary fortran preprocessing
                    if self.warp.cgnsgrid.cgnsstructured:
                        self.warp.processstructuredpatches()
                    else:
                        self.warp.processunstructuredpatches()

                    fullConn = self.warp.cgnsgrid.surfaceconn - 1
                    fullPts = self.warp.cgnsgrid.surfacepoints

                    nPatch = self.warp.cgnsgrid.getnpatch()
                    fullPatchNames = []
                    for i in range(nPatch):
                        fullPatchNames.append(self.warp.cgnsgrid.getsurf(i + 1).strip().lower())
                    symmFamilies = set()
                    if self.getOption("symmetrySurfaces") is None:
                        # Use all symmetry surfaces:
                        for i in range(len(self.warp.cgnsgrid.surfaceissymm)):
                            if self.warp.cgnsgrid.surfaceissymm[i]:
                                symmFamilies.add(fullPatchNames[i].lower())
                    else:
                        # The user has supplied a list of surface families
                        for name in self.getOption("symmetrySurfaces"):
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
                            curNodeSize = fullPatchSizes[i][0] * fullPatchSizes[i][1]
                            curCellSize = (fullPatchSizes[i][0] - 1) * (fullPatchSizes[i][1] - 1)

                            if fullPatchNames[i] in symmFamilies:
                                # Keep track of the families we've actually used
                                usedFams.add(fullPatchNames[i])

                                # Determine the average normal for these points:
                                conn = (
                                    fullConn[curCellIndex : curCellIndex + curCellSize * 4] - fullConn[curCellIndex] + 1
                                )
                                pts = fullPts[curNodeIndex : curNodeIndex + curNodeSize * 3]
                                avgNorm = self.warp.averagenormal(pts, conn, 4 * np.ones(curCellSize, "intc"))
                                planes.append([pts[0:3], avgNorm])
                            else:
                                # If we skipped, we increment the offset
                                curOffset += curNodeSize

                            curNodeIndex += curNodeSize * 3
                            curCellIndex += curCellSize * 4

                    # end for (root proc)

                else:  # unstructured

                    # We won't do this in general. The issue is that
                    # each of the elements needs to be checked
                    # individually since one sym BC may have multiple
                    # actual symmetry planes. This could be done, but
                    # since plane elemination code is in python and
                    # slow, we'll just defer this and make the user
                    # supply the symmetry planes.

                    raise Error(
                        "Automatic determine of symmetry surfaces is "
                        "not supported for unstructured CGNS. Please "
                        "specify the symmetry planes using the "
                        "'symmetryPlanes' option. See the _setSymmetryConditions()"
                        " documentation string for the option prototype."
                    )

                # Check if all supplied family names were actually
                # used. The user probably wants to know if a family
                # name was specified incorrectly.
                if self.comm.rank == 0:
                    if usedFams < symmFamilies:
                        missing = list(symmFamilies - usedFams)
                        Warning(
                            "Not all specified symm families that "
                            "were given were found the CGNS file. "
                            "The families not found are %s." % (repr(missing))
                        )

            elif self.fileType in ["OpenFOAM", "PLOT3D"]:

                # We could probably implement this at some point, but
                # it is not critical

                raise Error(
                    "Automatic determine of symmetry surfaces is "
                    "not supported for OpenFOAM or PLOT3D meshes. Please "
                    "specify the symmetry planes using the "
                    "'symmetryPlanes' option. See the _setSymmetryConditions()"
                    " documentation string for the option prototype."
                )

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
                d1 = p2 - np.dot(d, n1) * n1

                if np.linalg.norm(d1 - p2) / (np.linalg.norm(d) + 1e-30) > 1e-8:
                    return False

                return True

            uniquePlanes = []
            flagged = np.zeros(len(planes), "intc")
            for i in range(len(planes)):
                if not flagged[i]:
                    uniquePlanes.append(planes[i])
                    curPlane = planes[i]
                    # Loop over remainder to check:
                    for j in range(i + 1, len(planes)):
                        if checkPlane(curPlane[0], curPlane[1], planes[j][0], planes[j][1]):
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
                pstar = p2 - np.dot(d, n) * n

                normals.append(n)
                pts.append(pstar)

        normals = self.comm.bcast(np.array(normals))
        pts = self.comm.bcast(np.array(pts))

        # Let the user know what they are:
        if self.comm.rank == 0:
            print("+-------------------- Symmetry Planes -------------------+")
            print("|           Point                        Normal          |")
            for i in range(len(pts)):
                print(
                    "| (%7.3f %7.3f %7.3f)    (%7.3f %7.3f %7.3f) |"
                    % (
                        np.real(pts[i][0]),
                        np.real(pts[i][1]),
                        np.real(pts[i][2]),
                        np.real(normals[i][0]),
                        np.real(normals[i][1]),
                        np.real(normals[i][2]),
                    )
                )
            print("+--------------------------------------------------------+")

        # Now set the data into fortran.
        self.warp.setsymmetryplanes(pts.T, normals.T)

    # =========================================================================
    #                  Local OpenFOAM Functions
    # =========================================================================

    def _computeOFConn(self):

        """
        The user has specified an OpenFOAM mesh. Loop through the mesh data and
        create the arrays necessary to initialize the warping.
        """

        # Finally create the underlying data structure:
        faceSizes = []
        faceConn = []
        pts = []
        x0 = self.OFData["x0"]
        faces = self.OFData["faces"]
        connCount = 0

        for bName in self.OFData["boundaries"]:
            bType = self.OFData["boundaries"][bName]["type"].lower()

            if bType in ["patch", "wall", "slip"]:
                # Apparently OpenFOAM will list boundaries with zero
                # faces on them:
                nFace = len(self.OFData["boundaries"][bName]["faces"])
                if nFace > 0:
                    for iFace in self.OFData["boundaries"][bName]["faces"]:
                        face = faces[iFace]
                        nNodes = len(face)
                        # Pull out the 'len(face)' nodes from x0.
                        for j in range(nNodes):
                            pts.append(x0[face[j]])
                        faceSizes.append(nNodes)
                        faceConn.extend(range(connCount, connCount + nNodes))
                        connCount += nNodes

        pts = np.array(pts)
        faceConn = np.array(faceConn)
        return faceSizes, faceConn, pts

    def _readOFGrid(self, caseDir):
        """
        Read in the mesh points and connectivity from the polyMesh
        directory.

        NOTE: To ensure a consistant starting point in DAFoam. One needs to copy
        points_orig to points in pyDAFoam before initializing USMesh. This function
        will no longer perform this copying!

        Parameters
        ----------
        caseDir : str
            The directory containing the OpenFOAM mesh files
        """

        # import the pyOFM repo
        from pyofm import PYOFM

        ofm = PYOFM(comm=self.comm)

        # Copy the reference points file to points to ensure
        # consistant starting point
        self.OFData = ofm.getFileNames(caseDir, comm=self.comm)

        # Read in the volume points
        self.OFData["x0"] = ofm.readVolumeMeshPoints()

        # Read the face info for the mesh
        self.OFData["faces"] = ofm.readFaceInfo()

        # Read the boundary info
        self.OFData["boundaries"] = ofm.readBoundaryInfo(self.OFData["faces"])

        # Read the cell info for the mesh
        self.OFData["owner"], self.OFData["neighbour"] = ofm.readCellInfo()

        # Create the internal structures for volume mesh
        x0 = self.OFData["x0"]

        self.warp.createcommongrid(x0.T)

        ofm = None

    # =========================================================================
    #                     Utility functions
    # =========================================================================

    def _setMeshOptions(self):
        """Private function to set the options currently in
        self.options to the corresponding place in Fortran"""

        self.warp.gridinput.ldeffact = self.getOption("LdefFact")
        self.warp.gridinput.alpha = self.getOption("alpha")
        self.warp.gridinput.aexp = self.getOption("aExp")
        self.warp.gridinput.bexp = self.getOption("bExp")
        self.warp.gridinput.symmtol = self.getOption("symmTol")
        self.warp.gridinput.userotations = self.getOption("useRotations")
        self.warp.gridinput.zerocornerrotations = self.getOption("zeroCornerRotations")
        self.warp.gridinput.cornerangle = self.getOption("cornerAngle")
        self.warp.gridinput.errtol = self.getOption("errTol")
        self.warp.kd_tree.bucket_size = self.getOption("bucketSize")
        if self.getOption("evalMode") == "fast":
            self.warp.gridinput.evalmode = self.warp.gridinput.eval_fast
        elif self.getOption("evalMode") == "exact":
            self.warp.gridinput.evalmode = self.warp.gridinput.eval_exact

    def _processFortranStringArray(self, strArray):
        """Getting arrays of strings out of Fortran can be kinda nasty. This
        takes the array and returns a nice python list of strings"""
        shp = strArray.shape
        arr = strArray.reshape((shp[1], shp[0]), order="F")
        tmp = []

        for i in range(arr.shape[1]):
            tmp.append("")
            for j in range(arr.shape[0]):
                tmp[-1] += arr[j, i].decode()
            tmp[-1] = tmp[-1].strip().lower()

        return tmp

    def __del__(self):
        """Release all the mesh warping memory. This should be called
        automatically when the object is garbage collected."""
        self.warp.releasememory()
