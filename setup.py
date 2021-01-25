from setuptools import setup
import re

__version__ = re.findall(
    r"""__version__ = ["']+([0-9\.]*)["']+""",
    open("idwarp/__init__.py").read(),
)[0]

setup(
    name="idwarp",
    version=__version__,
    description="idwarp is a package deforming volume meshes with derivatives for optimization",
    long_description="""

      # IDWarp
      [![Build Status](https://travis-ci.com/mdolab/idwarp.svg?branch=master)](https://travis-ci.com/mdolab/idwarp)
      [![Documentation Status](https://readthedocs.com/projects/mdolab-idwarp/badge/?version=latest)](https://mdolab-idwarp.readthedocs-hosted.com/en/latest/?badge=latest)


      IDWarp uses an inverse distance method to modify the location of mesh volume nodes given a perturbation of the surface nodes.

      ## Documentation

      Please see the [documentation](https://mdolab-idwarp.readthedocs-hosted.com/en/latest/) for installation details and API documentation.

      To locally build the documentation, enter the `doc` folder and enter `make html` in terminal.
      You can then view the built documentation in the `_build` folder.


      ## Citation

      IDWarp is based on the theory presented in [this journal article](https://doi.org/10.1016/j.jcp.2011.09.021).
      """,
    long_description_content_type="text/markdown",
    keywords="mesh-warping warping mesh mesh-deformation optimization",
    author="",
    author_email="",
    url="https://github.com/mdolab/idwarp",
    license="Apache License Version 2.0",
    packages=[
        "idwarp",
    ],
    package_data={"idwarp": ["*.so"]},
    install_requires=[
        "numpy>=1.16",
        "petsc4py>=3.11",
        "mpi4py>=3.0",
        "mdolab-baseclasses>=1.4",
    ],
    extras_require={"testing": ["parameterized", "testflo"]},
    classifiers=["Operating System :: Linux", "Programming Language :: Python, Fortran"],
)
