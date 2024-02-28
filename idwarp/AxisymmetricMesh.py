import os
import shutil
import warnings
from dataclasses import dataclass

import numpy as np
from baseclasses import BaseSolver
from baseclasses.utils import Error
from mpi4py import MPI
from numpy.typing import ArrayLike

from .MExt import MExt


@dataclass(eq=True, unsafe_hash=True, frozen=True)
class MeshPatch:
    idxStart: int
    idxEnd: int
    point: ArrayLike
    normal: ArrayLike
    nodeSize: int
    cellSize: int
    patchName: str


class AxisymmetricMesh(BaseSolver):
    def __init__(self, options, comm=None, debug=False):
        """
        Create the AxisymmetricMesh object.

        Parameters
        ----------
        options :  dictionary
            A dictionary containing the options for the mesh movement
            strategies.

        comm : MPI_INTRA_COMM
            MPI communication (as obtained from mpi4py) on which to create
            the AxisymmetricMesh object. If not provided, MPI_COMM_WORLD is used.

        debug : bool
            Flag specifying if the MExt import is automatically deleted.
            This needs to be true ONLY when a symbolic debugger is used.
        """
        name = "IDWarp"
        category = "Volume mesh warping"

        if comm is None:
            comm = MPI.COMM_WORLD

        # Get the default options
        defOpts = self._getDefaultOptions()

        # Initialize the base class
        super().__init__(name, category, defaultOptions=defOpts, options=options, comm=comm)

        self.printOptions()

        # Check if warp was already set (complex inheritance issue)
        if not hasattr(self, "warp"):
            curDir = os.path.basename(os.path.dirname(os.path.realpath(__file__)))
            self.warp = MExt("libidwarp", curDir, debug=debug)._module

        # Initialize PETSc
        self.warp.initpetsc(self.comm.py2f())

        # Set the dtype to 'd'
        self.dtype = "d"

        # Set the fortran option values
        self._setMeshOptions()

        # Initialize class attributes
        self.warpInitialized = False
        self.fileType = self.getOption("fileType")
        fileName = self.getOption("gridFile")

        if not self.fileType == "CGNS":
            raise Error("Only CGNS file type is supported for axisymmetric mesh warping")

        self.warp.readcgns(fileName)

        if not self.warp.cgnsgrid.cgnsstructured:
            raise Error("Only structured CGNS files are supported for axisymmetric mesh warping")

        # Check the boco types on the root proc
        if self.comm.rank == 0:
            if 2 not in self.warp.cgnsgrid.bocoTypes and 17 not in self.warp.cgnsgrid.bocoTypes:
                raise Error("No 'wedge' nor 'symmPolar' boundary condition found in the CGNS file")

    @staticmethod
    def _getDefaultOptions():
        defOpts = {
            # General options
            "gridFile": [(str, type(None)), None],
            "fileType": [str, ["CGNS"]],
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
            # Axisymmetric only options
            "axiSymmetric": [bool, False],
            "axiAngle": [float, 0.0],
            "axiAxis": [list, [1.0, 0.0, 0.0]],
            "axiPlane": [(list, type(None)), None],
            "axiMirrorFamily": [(str, type(None)), None],
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

        # Get the restart file option
        restartFile = "" if self.getOption("restartFile") is None else self.getOption("restartFile")

        self._setSymmetryConditionsAxisymm()

        allPts = np.vstack(self.comm.allgather(pts))
        allFaceSizes = np.hstack(self.comm.allgather(faceSizes))

        # Be careful with conn, need to increment by the offset from the points
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
            dXvWarp = np.random.random(self.warp.griddata.warpmeshdof)  # , dtype=self.dtype)
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
        if self.fileType == "CGNS":
            if fileName is None:
                raise Error("fileName is not optional for writeGrid with gridType of CGNS")

            # Copy the default and then write
            if self.comm.rank == 0:
                shutil.copy(self.getOption("gridFile"), fileName)
            self.warp.writecgns(fileName)

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
            warnings.warn(
                "Using internally generated IDWarp surfaces. If "
                "this mesh object is to be used with an "
                "external solver, ensure the mesh object is "
                "passed to the solver immediatedly after it is created. "
                "The external solver must then call "
                "'setExternalMeshIndices()' and 'setSurfaceDefinition()' "
                "routines.",
                stacklevel=2,
            )

        conn = []
        pts = []
        faceSizes = []

        if self.comm.rank == 0:
            # Do the necessary fortran preprocessing
            self.warp.processstructuredpatches()

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

            # Check if all supplied family names were actually
            # used. The user probably wants to know if a family
            # name was specified incorrectly.
            if self.comm.rank == 0:
                if usedFams < surfaceFamilies:
                    missing = list(surfaceFamilies - usedFams)
                    warnings.warn(
                        "Not all specified surface families that "
                        "were given were found the CGNS file. "
                        "The families not found are %s." % (repr(missing)),
                        stacklevel=2,
                    )
                if len(usedFams) == 0:
                    raise Error(
                        "No surfaces were selected. Check the names "
                        "given in the 'specifiedSurface' option. The "
                        "list of families is %s." % (repr(list(fullPatchNames)))
                    )

    def _getFullPatchData(self):
        self.warp.processstructuredpatches()

        fullConn = self.warp.cgnsgrid.surfaceconn - 1
        fullPts = self.warp.cgnsgrid.surfacepoints

        nPatch = self.warp.cgnsgrid.getnpatch()
        fullPatchNames = [self.warp.cgnsgrid.getsurf(i + 1).strip().lower() for i in range(nPatch)]
        fullPatchNames = [val.decode("utf-8") if isinstance(val, bytes) else val for val in fullPatchNames]

        return fullPatchNames, fullConn, fullPts

    def _getSymmetryFamilies(self, fullPatchNames):
        symmFamilies = set()
        symmFamilies.update([fullPatchNames[i] for i, isSymm in enumerate(self.warp.cgnsgrid.surfaceissymm) if isSymm])

        # Decode family names that are bytes to convert them to strings
        symmFamilies = {val.decode("utf-8") if isinstance(val, bytes) else val for val in symmFamilies}

        return symmFamilies

    def _getSymmetryPlanesFromGeometry(self):
        fullPatchNames, fullConn, fullPts = self._getFullPatchData()
        symmFamilies = self._getSymmetryFamilies(fullPatchNames)
        symmPatches = []
        usedFams = set()

        fullPatchSizes = self.warp.cgnsgrid.surfacesizes.T

        curNodeIndex = 0
        curCellIndex = 0
        curOffset = 0

        for i, patchName in enumerate(fullPatchNames):
            # Get the node and cell sizes
            curNodeSize = fullPatchSizes[i][0] * fullPatchSizes[i][1]
            curCellSize = (fullPatchSizes[i][0] - 1) * (fullPatchSizes[i][1] - 1)

            if patchName in symmFamilies:
                # Keep track of the families we've actually used
                usedFams.add(patchName)

                # Determine the average normal for this patch
                conn1 = fullConn[curCellIndex : curCellIndex + curCellSize * 4]
                conn2 = fullConn[curCellIndex]
                conn = conn1 - conn2 + 1

                idxStart = curNodeIndex
                idxEnd = idxStart + curNodeSize * 3
                pts = fullPts[curNodeIndex : curNodeIndex + curNodeSize * 3]
                avgNorm = self.warp.averagenormal(pts, conn, 4 * np.ones(curCellSize, "intc"))

                # Decode patchName if it's a byte to convert it to a string
                patchName = patchName.decode("utf-8") if isinstance(patchName, bytes) else patchName
                patch = MeshPatch(idxStart, idxEnd, pts[:3], avgNorm, curNodeSize, curCellSize, patchName)
                symmPatches.append(patch)
            else:
                curOffset += curNodeSize

            curNodeIndex += curNodeSize * 3
            curCellIndex += curCellSize * 4

        self._checkUsedFamilies(usedFams, symmFamilies)

        return symmPatches, usedFams

    def _checkUsedFamilies(self, usedFams, symmFamilies):
        if usedFams < symmFamilies:
            missing = list(symmFamilies - usedFams)
            warnings.warn(
                "Not all specified symm families that "
                "were given were found the CGNS file. "
                "The families not found are %s." % (repr(missing)),
                stacklevel=2,
            )

    def _printSymmetryPlanes(self, pts, normals):
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

    def _setSymmetryConditionsAxisymm(self):
        # Initialize these variables on all procs
        pts = []
        normals = []
        mirrorFamName = None
        if self.comm.rank == 0:
            mirrorFamName = self.getOption("axiMirrorFamily")

            if mirrorFamName is None:
                raise Error("Must specify an 'axiMirrorFamily' for axisymmetric mesh warping.")

            # We need the information for all symm planes for axisymmetric warping
            symmPatches, usedFams = self._getSymmetryPlanesFromGeometry()

            # Should only have two symmetry planes for an axisymm mesh
            if len(usedFams) != 2:
                raise Error("Could not detect 2 unique symmetry planes for an axisymmetric mesh.")

            mirrorPatches = [patch for patch in symmPatches if patch.patchName == mirrorFamName]
            rotPatches = [patch for patch in symmPatches if patch.patchName != mirrorFamName]

            # Get the first two patches to represent the symmetry planes
            plane1 = mirrorPatches[0]
            plane2 = rotPatches[0]

            pts.append(plane1.point)
            normals.append(plane1.normal)

            self._printSymmetryPlanes([plane1.point, plane2.point], [plane1.normal, plane2.normal])

        # Broadcast the mirror family name to all procs
        mirrorFamName = self.comm.bcast(mirrorFamName)
        self.warp.griddata.setmirrorfamily(mirrorFamName)

        # We only want to set the mirror plane as the symmetry plane in the
        # fortran layer.  This "tricks" the mesh warping algorithm so that
        # the mirroring will get equal weights to preserve the symmetry plane
        # warping.
        pts = self.comm.bcast(np.array(pts))
        normals = self.comm.bcast(np.array(normals))
        axiPlane = self.getOption("axiPlane")
        if axiPlane is None:
            raise Error(
                "Must supply a point and normal vector using the "
                "'axiPlane' option to perform axisymmetric mesh warping."
            )

        # Set a "phantom" plane that preserves zero-warping along the degenerate line BC
        pts = np.vstack((pts, axiPlane[0]))
        normals = np.vstack((normals, axiPlane[1]))

        self.warp.setsymmetryplanes(pts.T, normals.T)

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
        self.warp.griddata.axisymm = self.getOption("axiSymmetric")
        if self.getOption("evalMode") == "fast":
            self.warp.gridinput.evalmode = self.warp.gridinput.eval_fast
        elif self.getOption("evalMode") == "exact":
            self.warp.gridinput.evalmode = self.warp.gridinput.eval_exact

        if self.getOption("axiSymmetric"):
            if self.getOption("axiAngle") != 0.0:
                angle = self.getOption("axiAngle")
                self.warp.griddata.axisymmangle = np.deg2rad(angle)  # Convert to radians
            else:
                raise Error("Must supply an 'axiSymmAngle' for axisymmetric mesh warping.")

            if any([i != 0.0 for i in self.getOption("axiAxis")]):
                self.warp.griddata.axisymmaxis = np.array(self.getOption("axiAxis"), self.dtype)
            else:
                raise Error("Must supply an 'axiAxis' for axisymmetric mesh warping.")

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
