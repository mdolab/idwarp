#!/usr/bin/python
"""
This is the complex version of the MultiUSMesh class.
See MultiUnstructuredMesh.py for more info.
"""

# =============================================================================
# Imports
# =============================================================================
import os
from .MultiUnstructuredMesh import MultiUSMesh
from . import MExt

# =============================================================================
# USMesh_C class
# =============================================================================


class MultiUSMesh_C(MultiUSMesh):
    """
    Represents a (complex) multi-component mesh
    """

    def __init__(self, *args, **kwargs):
        """Initialize the object."""

        debug = False
        if "debug" in kwargs:
            debug = kwargs["debug"]

        curDir = os.path.basename(os.path.dirname(os.path.realpath(__file__)))
        self.warp = MExt.MExt("libidwarp_cs", curDir, debug=debug)._module
        MultiUSMesh.__init__(self, *args, dtype="D", **kwargs)
