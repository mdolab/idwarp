from setuptools import setup
import re
import os

__version__ = re.findall(
    r"""__version__ = ["']+([0-9\.]*)["']+""",
    open("idwarp/__init__.py").read(),
)[0]

this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="idwarp",
    version=__version__,
    description="idwarp is a package deforming volume meshes with derivatives for optimization",
    long_description=long_description,
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
