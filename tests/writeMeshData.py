"""
This script was used to generate the onera_m6.json input file needed for the MultiUSMesh test.
It could also be useful for creating other tests in the future.
We use ADflow to get mesh information and then write it as a dictionary in a JSON file.
Unlike USMesh, MultiUSMesh needs surface information before generating internal surfaces.
By saving the mesh data, we avoid adding ADflow as a testing dependency.
"""

import os
from adflow import ADFLOW
from baseclasses.utils import writeJSON

baseDir = os.path.dirname(os.path.abspath(__file__))
gridFile = os.path.join(baseDir, "../input_files/onera_m6.cgns")

# Set up ADflow
aeroOptions = {
    "gridFile": gridFile,
    "outputDirectory": baseDir,
    "MGCycle": "sg",  # This is needed to avoid problems with coarse grid initialization
}
CFDSolver = ADFLOW(options=aeroOptions)

# Get the mesh data from ADflow
meshInd = CFDSolver.getSolverMeshIndices()
conn, faceSizes, cgnsBlockIDs = CFDSolver.getSurfaceConnectivity(
    CFDSolver.meshFamilyGroup, includeZipper=False, includeCGNS=True
)
pts = CFDSolver.getSurfaceCoordinates(CFDSolver.meshFamilyGroup, includeZipper=False)

# Save the ADflow output
meshData = {
    "meshInd": meshInd,
    "conn": conn,
    "faceSizes": faceSizes,
    "cgnsBlockIDs": cgnsBlockIDs,
    "pts": pts,
}
jsonFile = os.path.join(baseDir, "onera_m6.json")
writeJSON(jsonFile, meshData)
