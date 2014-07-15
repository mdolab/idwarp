MODULE_FILES="../modules/communication.f90 ../modules/precision.f90"

IO_FILES="../IO/readUnstructuredCGNS.f90"

GEOCALC_FILES="../geoCalcs/getSurfaceElementCenters.f90"

f2py -m warpustruct -h warpustruct.pyf $IO_FILES $MODULE_FILES $GEOCALC_FILES