.. _pyWarpUstruct_building:

Building
--------

All the core computations in :ref:`pyWarpUstruct` are coded in Fortran.  It
is therefore necessary to build this library before using
:ref:`pyWarpUstruct`.

To see a list of architectures that :ref:`pyWarpUstruct` has been known to
compile on run::
   
   make

from the root directory. 

The easiest approach to to try the closest one to your system and
attempt a build using (for example)::

   make LINUX_INTEL_OPENMPI

If everything was successful, the following lines will be printed to
the screen (near the end)::

   Testing if module warp can be imported...
   Module warp was successfully imported.

If you don't see this, it will be necessary to configure the build
manually. To configure manually, first copy a default configuration
file from the defaults folder like this (run this in the root
directory)::
  
   cp config/defaults/config.LINUX_INTEL_OPENMPI.mk config

Now open ``config/config.LINUX_INTEL_OPENMPI.mk`` which should look like::

  # ----------------------------------------------------------------------
  # Config file for Intel ifort  with OpenMPI
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
  FF90_GEN_FLAGS = 
  CC_GEN_FLAGS   =

  FF90_OPT_FLAGS   =  -fPIC -r8 -O2  
  CC_OPT_FLAGS     = -O2 -fPIC

  # ------- Define Linker Flags ------------------------------------------
  LINKER_FLAGS = -nofor_main

  # ------- Define Petsc Info --- Should not need to modify this -----
  include ${PETSC_DIR}/conf/variables
  PETSC_INCLUDE_FLAGS=${PETSC_CC_INCLUDES} -I$(PETSC_DIR)
  PETSC_LINKER_FLAGS=${PETSC_LIB}

It will most likely be necessary to modify the ``CGNS_INCLUDE_FLAGS``
and the ``CGNS_LINKER_FLAGS`` variables. It is also necessary to
``PETSc`` already compiled including support for
``SuperLU_dist``. After changes to the configuration file, run ``make
clean`` before attempting a new build. 
