# Config File for LINUX and INTEL Compiler
AR       = ar
AR_FLAGS = -rvs
RM       = /bin/rm -rf

# Fortran compiler and flags
FF90        = mpif90
FF90_FLAGS  = -fdefault-real-8 -O2 -fPIC

# C compiler and flags
CC       = mpicc
CC_FLAGS   = -O2 -fPIC

# Define potentially different python, python-config and f2py executables:
PYTHON = python
PYTHON-CONFIG = python-config
F2PY = f2py

# Define additional flags for linking
LINKER_FLAGS = -L /usr/lib -llapack

# ------- Define Petsc Info --- Should not need to modify this -----
include ${PETSC_DIR}/conf/variables
PETSC_INCLUDE_FLAGS=${PETSC_CC_INCLUDES} -I$(PETSC_DIR)
PETSC_LINKER_FLAGS=${PETSC_LIB}

# ------- Define CGNS Inlcude and linker flags -------------------------
CGNS_INCLUDE_FLAGS=-I$(HOME)/packages/cgnslib_3.2.1/src
CGNS_LINKER_FLAGS=-L$(HOME)/packages/cgnslib_3.2.1/src -lcgns