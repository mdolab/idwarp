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
        curDir = os.path.dirname(os.path.realpath(__file__))
        self.warp = MExt.MExt('warpustruct_cs', [curDir])._module
        USMesh.__init__(self, *args, **kwargs)
        self.dtype = 'D'
