# Config File for LINUX and INTEL Compiler
AR       = ar
AR_FLAGS = -rvs
RM       = /bin/rm -rf

# Fortran compiler and flags
FF90        = ifort
FF90_FLAGS  = -r8 -O2 -fPIC -DUSE_NO_PETSC

# C compiler and flags
CC       = gcc
CC_FLAGS   = -O2 -fPIC

# Define potentially different python, python-config and f2py executables:
PYTHON = python
PYTHON-CONFIG = python-config
F2PY = f2py

# Define additional flags for linking
LINKER_FLAGS = -L /usr/lib -llapack

# ------- Define CGNS Inlcude and linker flags -------------------------
CGNS_INCLUDE_FLAGS=-I$(HOME)/packages/cgnslib_3.2.1/src
CGNS_LINKER_FLAGS=-L$(HOME)/packages/cgnslib_3.2.1/src -lcgns