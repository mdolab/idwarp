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

To see a list of architectures that IDWarp has been known to compile on run

.. prompt:: bash

   make

from the root directory.

Follow the instructions and copy the closest architecture file and
attempt a build using

.. prompt:: bash

   make

If everything was successful, the following lines will be printed to
the screen (near the end)::

   Testing if module idwarp can be imported...
   Module idwarp was successfully imported.

If you don't see this, it will be necessary modify the configure
options in the config file.

It will most likely be necessary to modify the ``CGNS_INCLUDE_FLAGS``
and the ``CGNS_LINKER_FLAGS`` variables. After changes to the
configuration file, run ``make clean`` before attempting a new build.

Lastly, to build the Python interface, go to the root directory and type

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

Once the library is built run the following to install the python module and library into your environment.

.. prompt:: bash

   pip install .[complex]