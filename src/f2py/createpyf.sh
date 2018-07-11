#MODULE_FILES="../modules/communication.f90  ../modules/gridData.F90 ../modules/gridInput.f90 ../modules/kd_tree.F90 ../modules/pointReduce.F90 ../modules/structuredCGNS.F90"
MODULE_FILES=" ../modules/gridData.F90 ../modules/gridInput.f90" # ../modules/kd_tree.F90"

#IO_FILES="../IO/createGrid.F90 ../IO/patchIO.f90 ../IO/readStructuredCGNS.F90 ../IO/readUnstructuredCGNS.f90 ../IO/writeStructuredCGNS.F90 ../IO/writeUnstructuredCGNS.f90"
IO_FILES="../IO/createGrid.F90 ../IO/patchIO.F90 ../IO/readStructuredCGNS.F90  ../IO/readUnstructuredCGNS.f90 ../IO/writeStructuredCGNS.F90 ../IO/writeUnstructuredCGNS.f90"

WARP_FILES="../warp/getCommonVolumeCoordinates.F90 ../warp/getdXs.F90 ../warp/getSurfaceCoordinates.F90 ../warp/getVolumeCoordinates.F90 ../warp/initializeWarping.F90 ../warp/setdXs.F90 ../warp/setExternalMeshIndices.F90 ../warp/setSurfaceCoordinates.F90 ../warp/warpDeriv.F90 ../warp/warpMesh.F90"

f2py -m idwarp -h test.pyf $IO_FILES $MODULE_FILES  $WARP_FILES
