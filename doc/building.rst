.. _IDWarp_building:

Building
--------
IDWarp depends on the follow libraries:
- CGNS Library
- PETSc
- MPI

See the MDO Lab installation guide :ref:`here <mach-aero:installThirdPartyPackages>` for the supported versions and installation instructions.

All the core computations in IDWarp are coded in Fortran.
It is therefore necessary to build this library before using IDWarp.

To see a list of architectures that IDWarp has been known to
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

It will most likely be necessary to modify the ``CGNS_INCLUDE_FLAGS``
and the ``CGNS_LINKER_FLAGS`` variables. After changes to the
configuration file, run ``make clean`` before attempting a new build.

Lastly, to build the Python interface, go to the root directory and type::

   pip install .
