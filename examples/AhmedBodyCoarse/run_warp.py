# a basic script to run the case in this directory
import os
from mpi4py import MPI
from pygeo import DVGeometry
from pyspline import Curve
from idwarp import USMesh
import numpy

gcomm = MPI.COMM_WORLD

meshOptions = {
    "gridFile": os.getcwd(),
    "fileType": "openFoam",
    "symmetryPlanes": [[[0, 0, 0], [0, 1, 0]]],
    "aExp": 3,
    "bExp": 5,
    "alpha": 1.0,
    "LdefFact": 0.20,
}

mesh = USMesh(options=meshOptions, comm=gcomm)

coords0 = mesh.getSurfaceCoordinates()

# setup FFD
FFDFile = "./FFD/globalFFD.fmt"
DVGeo = DVGeometry(FFDFile)
# Setup curves for ref_axis
x = [-2.0, 0.0, 0.1, 1.044, 5.0]
y = [0.1, 0.1, 0.1, 0.1, 0.1]
z = [0.1, 0.1, 0.1, 0.1, 0.1]

nLength = len(x)

c1 = Curve(x=x, y=y, z=z, k=2)
DVGeo.addRefAxis("bodyAxis", curve=c1, axis="z")

DVGeoChild = DVGeometry("./FFD/bodyFittedFFD.fmt", child=True)

# Setup curves for ref_axis
x1 = [0.0, 0.1, 0.862, 1.044]
y1 = [0.1, 0.1, 0.1, 0.1]
z1 = [0.194, 0.194, 0.194, 0.13]
# z1 = [0.338,0.338,0.338,0.21]
# z1 = [0.338,0.338,0.338,0.338]

nLengthChild = len(x1)

c2 = Curve(x=x1, y=y1, z=z1, k=2)
DVGeoChild.addRefAxis("localBodyAxis", curve=c2, axis="z")


def rampAngle(val, geo):
    C = geo.extractCoef("localBodyAxis")

    # the value will be ramp angle in degree.
    # start with a conversion to rads
    angle = (val[0]) * numpy.pi / 180.0

    # Overall length needs to stay a 1.044, so use that as a ref for
    # the final mesh point

    # set the target length
    lTarget = 0.222
    hInit = 0.21 - 0.05

    # compute the coefficient deltas
    dx = lTarget * numpy.cos(angle)
    dz = lTarget * numpy.sin(angle)

    topEdge = 0.338 - dz
    rearHeight = topEdge - 0.05
    coefPoint = rearHeight / 2.0 + 0.05
    scalez = rearHeight / hInit

    # Set the coefficients
    C[3, 0] = 1.044
    C[2, 0] = C[3, 0] - dx
    C[2, 2] = 0.194
    C[3, 2] = coefPoint

    geo.restoreCoef(C, "localBodyAxis")

    geo.scale_z["localBodyAxis"].coef[3] = scalez

    return


def doubleRampAngle(val, geo):
    C = geo.extractCoef("localBodyAxis")

    # the values will be the upper and lower ramp angle
    # in degree. Start with a conversion to rads.
    upperAngle = (val[0]) * numpy.pi / 180.0
    lowerAngle = (val[1]) * numpy.pi / 180.0

    # Overall length needs to stay a 1.044, so use that as a ref for
    # the final mesh point

    # set the target length
    lTarget = 0.222
    hInit = 0.21 - 0.05

    # compute the coefficient deltas
    dx = lTarget * numpy.cos(upperAngle)
    dzUpper = lTarget * numpy.sin(upperAngle)
    dzLower = lTarget * numpy.sin(lowerAngle)

    topEdge = 0.338 - dzUpper
    rearHeight = topEdge - 0.05 + dzLower
    coefPoint = rearHeight / 2.0 + 0.05 - dzLower
    scalez = numpy.sqrt((rearHeight / hInit) ** 2)

    # Set the coefficients
    C[3, 0] = 1.044
    C[2, 0] = C[3, 0] - dx
    C[2, 2] = 0.194
    C[3, 2] = coefPoint

    geo.restoreCoef(C, "localBodyAxis")

    geo.scale_z["localBodyAxis"].coef[3] = scalez

    return


def length(val, geo):
    C = geo.extractCoef("bodyAxis")

    for i in range(len(C)):
        # print 'C',i,C[i,0],val[i]
        C[i, 0] = val[i]

    # end
    geo.restoreCoef(C, "bodyAxis")
    return


def angleVars(val, geo):
    C = geo.extractCoef("localBodyAxis")

    for i in range(len(C)):
        C[i, 2] = val[i]
    # end
    geo.restoreCoef(C, "localBodyAxis")


def noseLength(val, geo):
    C = geo.extractCoef("bodyAxis")

    length = val[0]

    C[1, 0] = C[2, 0] - length

    geo.restoreCoef(C, "bodyAxis")

    return


# lower = [-2.,-2.,-2.,-2.,5.]
# upper = [-2.,5.,5.,5.,5.]
# DVGeo.addGeoDVGlobal('length', x, length,
#                      lower=lower, upper=upper, scale=1.0)
DVGeo.addGeoDVGlobal("noseLength", 0.3, noseLength, lower=0, upper=0.5, scale=1.0)

# DVGeoChild.addGeoDVGlobal('rampAngle', 35.1, rampAngle,
#                      lower=0., upper=90., scale=1.0)
DVGeoChild.addGeoDVGlobal(
    "doubleRampAngle", [35.1, -5.0], doubleRampAngle, lower=[0.0, -45.0], upper=[45.0, 5.0], scale=[1.0, 1.0]
)
# lowerA = [0.,0.,0.,0.]
# upperA = [0.3,0.3,0.3,0.3]
# DVGeoChild.addGeoDVGlobal('angleVars', z1, angleVars,
#                      lower=lowerA, upper=upperA, scale=1.0)

# lowerL = [-1.,-1.,-1.,-1.]
# upperL = [2.0,2.0,2.0,2.0]
# DVGeoChild.addGeoDVGlobal('noseLen', x1, noseLength,
#                      lower=lowerL, upper=upperL, scale=1.0)

# DVGeo.addGeoDVGlobal('angleVars', z1, angleVars,
#                      lower=lowerA, upper=upperA, scale=1.0)

# Add the child to the parent
DVGeo.addChild(DVGeoChild)

ptSetName = "allSurfs"
freezeDict = {}  # '0':['jLow'],'1':['jLow'],'2':['jLow']}#'0':['jLow'],'1':['jHigh','jLow']}
DVGeo.addPointSet(coords0, ptSetName, faceFreeze=freezeDict)
xDV = DVGeo.getValues()

# Required design variables
# Rear ramp angle, fixed 200 mm length
# overall length
# nose length
# Ramp shape
# Ground separation
# Lower ramp angle

# Case 1: Rear Ramp angle, fixed length
# Case 2: Upper and lower ramp angles, fixed length
# Case 3: Nose length ( Do with global FFD)
# Case 4: Overall length ( Do with global FFD)
# Case 5: Shape

# xDV['length'][2] = 1.75#2.0#1.05
# xDV['angleVars'][2] = 0.15
# xDV['angleVars'][0] = 0.19
# xDV['noseLen'][0] = -0.1
# xDV['angleVars'][1] = 0.18
# xDV['angleVars'][2] = 0.18
# xDV['angleVars'][3] = 0.12

DVGeo.setDesignVars(xDV)
mesh.setSurfaceCoordinates(DVGeo.update(ptSetName))
mesh.warpMesh()
DVGeo.writeTecplot("warpedFFD.dat")
# mesh.writeOFGridTecplot('warped.dat')
mesh.writeGrid()

# # Repeat ================
# #xDV['length'][2] = 1.25#2.0#1.05
# xDV['angleVars'][2] = 0.3

# DVGeo.setDesignVars(xDV)
# #coords = DVGeo.update(ptSetName)
# # for i in range(coords0.shape[0]):
# #     if coords0[i,1]==0:
# #         print 'x',coords[i,:]
# #     # end
# # # end
# # for i in range(coords.shape[0]):
# #     print 'x',coords0[i,:],coords[i,:]
# # # print DVGeo.update(ptSetName)
# #sys.exit(0)
# mesh.setSurfaceCoordinates(DVGeo.update(ptSetName))
# DVGeo.writeTecplot('warpedFFD2.dat')

# mesh.warpMesh()
# #print 'mesh warped'
# mesh.writeOpenFOAMVolumePoints()
# #print 'points updated'
# meshName = os.path.join(os.getcwd(),"testAhmedMesh2")
# mesh.writeGridTecplot(meshName)
# #print 'file written'
