# ----------------------------------------------------------------------
# Config file for Gfortran  with OpenMPI
# ----------------------------------------------------------------------

# ------- Define a possible parallel make ------------------------------
PMAKE = make -j 4

# ------- Define the MPI Compilers--------------------------------------
FF90 = mpif90
CC   = mpicc

# ------- Define CGNS Inlcude and linker flags -------------------------
CGNS_INCLUDE_FLAGS = -I/usr/local/include
CGNS_LINKER_FLAGS = -Wl,-rpath,/usr/local/lib -lcgns

# ------- Define Compiler Flags ----------------------------------------
FF90_GEN_FLAGS = -fPIC
CC_GEN_FLAGS   = -fPIC

FF90_OPT_FLAGS   =  -fPIC -fdefault-real-8 -O2
CC_OPT_FLAGS     = -O2

# ------- Define Linker Flags ------------------------------------------
LINKER_FLAGS = 

# ------- Define Petsc Info --- Should not need to modify this -----
include ${PETSC_DIR}/conf/variables
PETSC_INCLUDE_FLAGS=${PETSC_CC_INCLUDES} -I$(PETSC_DIR)
PETSC_LINKER_FLAGS=${PETSC_LIB}

# Define potentially different python, python-config and f2py executables:
PYTHON = python
PYTHON-CONFIG = python-config
F2PY = f2py

