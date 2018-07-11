#!/usr/bin/python
"""
This is the complex version of the Unstructured mesh warping. See
UnstructuredMesh.py for more info.
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
    Represents a (Complex) Unstructured mesh
    """
    def __init__(self, *args, **kwargs):
        """Initialize the object."""

        debug = False
        if 'debug' in kwargs:
            debug = kwargs['debug']

        curDir = os.path.dirname(os.path.realpath(__file__))
        self.warp = MExt.MExt('warpustruct_cs', [curDir], debug=debug)._module
        MultiUSMesh.__init__(self, dtype='D', *args, **kwargs)
