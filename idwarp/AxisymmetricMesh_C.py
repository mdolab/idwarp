import os

from . import MExt
from .AxisymmetricMesh import AxisymmetricMesh


class AxisymmetricMesh_C(AxisymmetricMesh):
    """
    Represents a (Complex) Axisymmetric mesh
    """

    def __init__(self, *args, **kwargs):
        """Initialize the object."""

        debug = False
        if "debug" in kwargs:
            debug = kwargs["debug"]

        curDir = os.path.basename(os.path.dirname(os.path.realpath(__file__)))
        self.warp = MExt.MExt("libidwarp_cs", curDir, debug=debug)._module
        AxisymmetricMesh.__init__(self, *args, **kwargs)
        self.dtype = "D"
