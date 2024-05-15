# ----------------------------------------------------------------------
# Config file for Intel ifort
# ----------------------------------------------------------------------

# ------- Define a possible parallel make ------------------------------
PMAKE = make -j 4

# ------- Define the MPI Compilers--------------------------------------
ifdef I_MPI_ROOT # Using Intel MPI
  ifneq ($(shell command -v icc 2> /dev/null;),)
    # icc only exists on older Intel versions
    # Assume that we want to use the old compilers
    FF90 = mpiifort
    CC = mpiicc
  else
    # Use the new compilers
    FF90 = mpiifx
    CC = mpiicx
  endif
else # Using HPE MPI
  FF90 = ifort -lmpi
  CC   = icc -lmpi
endif

# ------- Define CGNS Inlcude and linker flags -------------------------
# Define the CGNS include directory and linking flags for the CGNS library.
# We are assuming that HDF5 came from PETSc so it is included in ${PETSC_LIB}.
# Otherwise you will have to specify the HDF5 library.
CGNS_INCLUDE_FLAGS=-I$(CGNS_HOME)/include
CGNS_LINKER_FLAGS=-L$(CGNS_HOME)/lib -lcgns

# ------- Define complexify inlcude and linker flags -------------------------
COMPLEXIFY_INCLUDE_FLAGS=-I$(COMPLEXIFY_DIR)/include
COMPLEXIFY_LINKER_FLAGS=-L$(COMPLEXIFY_DIR)/lib -lcomplexify

# ------- Define Compiler Flags ----------------------------------------
FF77_FLAGS = -fPIC -r8 -O2
FF90_FLAGS = ${FF77_FLAGS} -stand f08
C_FLAGS    = -fPIC -O2

# ------- Define Linker Flags ------------------------------------------
LINKER_FLAGS = -fPIC -nofor-main

# ------- Define Petsc Info --- Should not need to modify this -----
include ${PETSC_DIR}/lib/petsc/conf/variables
PETSC_INCLUDE_FLAGS=${PETSC_CC_INCLUDES} -I$(PETSC_DIR)
PETSC_LINKER_FLAGS=${PETSC_LIB}

# Define potentially different python, python-config and f2py executables:
PYTHON = python
PYTHON-CONFIG = python3-config # use python-config for python 2
F2PY = f2py
