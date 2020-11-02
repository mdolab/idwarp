#!/usr/bin/python
"""
This is the complex version of the Unstructured mesh warping. See
UnstructuredMesh.py for more info.
"""

# =============================================================================
# Imports
# =============================================================================
import os
from .UnstructuredMesh import USMesh
from . import MExt

# =============================================================================
# USMesh_C class
# =============================================================================


class USMesh_C(USMesh):
    """
    Represents a (Complex) Unstructured mesh
    """

    def __init__(self, *args, **kwargs):
        """Initialize the object."""

        debug = False
        if "debug" in kwargs:
            debug = kwargs["debug"]

        curDir = os.path.dirname(os.path.realpath(__file__))
        self.warp = MExt.MExt("idwarp_cs", [curDir], debug=debug)._module
        USMesh.__init__(self, *args, **kwargs)
        self.dtype = "D"
