.. _IDWarp_building:

Building
--------

All the core computations in :ref:`IDWarp` are coded in Fortran.  It
is therefore necessary to build this library before using
:ref:`IDWarp`.

To see a list of architectures that :ref:`IDWarp` has been known to
compile on run::
   
   make

from the root directory. 

Follow the instructions and copy the closest architecture file and
attempt a build using ::

   make

If everything was successful, the following lines will be printed to
the screen (near the end)::

   Testing if module idwarp can be imported...
   Module idwarp was successfully imported.

If you don't see this, it will be necessary modify the configure
options in the config file. 

Now open ``config/config.mk`` which should look like::

  # ----------------------------------------------------------------------
  # Config file for Gfortran  with OpenMPI
  # ----------------------------------------------------------------------

  # ------- Define a possible parallel make ------------------------------
  PMAKE = make -j 4

  # ------- Define the MPI Compilers--------------------------------------
  FF90 = mpif90
  CC   = mpicc

  # ------- Define CGNS Inlcude and linker flags -------------------------
  CGNS_INCLUDE_FLAGS=-I$(HOME)/packages/cgnslib_3.2.1/src
  CGNS_LINKER_FLAGS=-L$(HOME)/packages/cgnslib_3.2.1/src -lcgns

  # ------- Define Compiler Flags ----------------------------------------
  FF90_GEN_FLAGS = -fPIC
  CC_GEN_FLAGS   = -fPIC

  FF90_OPT_FLAGS   =  -fPIC -fdefault-real-8 -O2
  CC_OPT_FLAGS     = -O2

  # ------- Define Linker Flags ------------------------------------------
  LINKER_FLAGS = 

  # ------- Define Petsc Info --- Should not need to modify this -----
  include ${PETSC_DIR}/lib/petsc/conf/variables # PETSc 3.6+
  #include ${PETSC_DIR}/conf/variables # PETSc 3.5
  PETSC_INCLUDE_FLAGS=${PETSC_CC_INCLUDES} -I$(PETSC_DIR)
  PETSC_LINKER_FLAGS=${PETSC_LIB}

  # Define potentially different python, python-config and f2py executables:
  PYTHON = python
  PYTHON-CONFIG = python-config
  F2PY = f2py

It will most likely be necessary to modify the ``CGNS_INCLUDE_FLAGS``
and the ``CGNS_LINKER_FLAGS`` variables. After changes to the
configuration file, run ``make clean`` before attempting a new build.
