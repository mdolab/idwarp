.. _IDWarp_building:

Building
========

Requirements
------------
IDWarp depends on the follow libraries:

- CGNS Library
- PETSc
- MPI
- Complexify (see the `Complex build`_ section)

Building
--------

See the MDO Lab installation guide :ref:`here <mach-aero:installThirdPartyPackages>` for the supported versions and installation instructions.

All the core computations in IDWarp are coded in Fortran.
It is therefore necessary to build this library before using IDWarp.

To start, find a configuration file close to your current setup in::

    config/defaults

and copy it to ``config/config.mk``. For example:

.. prompt:: bash

    cp config/defaults/config.LINUX_GFORTRAN_OPENMPI.mk config/config.mk

Once you have copied the config file, compile ADflow by running:

.. prompt:: bash

   make

If everything was successful, the following lines will be printed to the screen (near the end)::

   Testing if module idwarp can be imported...
   Module idwarp was successfully imported.

If you don't see this, it will be necessary to configure the build manually.
To configure manually, open ``config/config.mk`` and modify options as necessary.
Common issues are often that dependency variable paths, such as ``CGNS_INCLUDE_FLAGS`` and ``CGNS_LINKER_FLAGS``, point to an incorrect location and need to be updated.
After changes to the configuration file, run ``make clean`` before attempting a new build.

Lastly, to build and install the Python interface, type

.. prompt:: bash

   pip install .

Some features require additional dependencies.
Using IDWarp with OpenFOAM meshes in DAFoam requires ``pyOFM``.
This dependency can be checked with

.. prompt:: bash

   pip install .[dafoam]

Using MultiUSMesh requires ``cgnsUtilities``, which can be checked with

.. prompt:: bash

   pip install .[multi]


Complex build
-------------
Its possible to build a "complexified" version of IDWarp directly from the real version.
To build IDWarp in complex mode, a complex version of PETSc and the `Complexify <https://github.com/mdolab/complexify>`__ module and library are required.
Once installed and configured, run the following to build the complexified library

.. prompt:: bash

   export PETSC_ARCH=$PETSC_ARCH_COMPLEX
   make -f Makefile_CS

Note that, ``PETSC_ARCH``, **must** be set and point to the complex PETSc before the code is compiled.
In the above example, an intermediate convenience variable, ``PETSC_ARCH_COMPLEX``, defines the complex PETSc arch path.
Once the library is built run the following to install the python module and library into your environment.

.. prompt:: bash

   pip install .[complex]