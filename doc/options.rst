.. _idwarp_options:

Options
=======

Here are the options currently available in IDWarp.

======================================  ==========  ===========================================   ================================================================================================================================================================================
Parameter                                  Type       Default                                       Description
======================================  ==========  ===========================================   ================================================================================================================================================================================
`gridFile`                              `str`       `None`                                        This is the grid file to use. It must always be spcified. It may be a structured or
                                                                                                  unstructured CGNS file or it may be an OpenFOAM directory containing a mesh specification.

`fileType`                              `str`       `CGNS`                                        Specify the type of grid file. Valid options are 'CGNS' or 'OpenFOAM'

`specifiedSurfaces`                     `list`      `None`                                        This option is used to specify which surfaces are used to build the surface tree where 
                                                                                                  deformations are to be specified. The default is 'None' which will automatically use all wall-type surfaces
                                                                                                  in the grid file. For CGNS files this corresponds to the following boundary condtiions:
                                                                                                  'BCWall', 'BCWallViscous', 'BCWallViscousHeatFlux', 'BCWallViscousAdiabatic', 'BCWallInviscid'.
                                                                                                  For OpenFOAM files, all 'Patch' and 'Wall' surfaces are assumed by default. If a non-None value
												  is given it should be list of families the use wants to use to generate the surface definition. 

`symmetrySurfaces`                      `list`      `None`                                        This option is used to specify which surfaces are used to determine symmetry planes. The default
                                                                                                  is 'None' which will automatically use the 'BCSymmetryPlane' contidions found the CGNS files. 
                                                                                                  This option is only valid for structured CGNS files. If a non-None value is given it should be 
												  a list of families. 

`symmetryPlanes`                        `list`      `None`                                        This is sort of a "last-resort" option. It is used to overwrite and explicitly define symmetry conditions
                                                                                                  IDWarp is to use. It is the only method for specifying symmetry for unstructured CGNS files and OpenFOAM files. 
												  For the 'symmetrySurfaces' option to be active, 'symmetryPlanes' must be None. If 'symmetryPlanes' is not 'None' it
												  is expected to be a list of the following form: [ [[x1, y1, z1], [nx1, ny1, nz1]], [[x2, y2, z2], [nx2, ny2, nz2]] ]. 
												  The previous example defines two symmetry planes using a point-normal approach. The first plane is defined by pt=(x1,y1,z1) with 
												  normal=(nx1, ny1, nz1) and the second plane is defined with pt=(x2,y2,z2), normal=(nx2, ny2, nz2). The normal direction should be 
												  normalized to unit magnitude. Note that no symmetry may be specified with the option 'symmetryPlanes':[]. 

`aExp`                                  `float`     `3.0`                                         This is the first exponent in the inverse distance calc. However, for efficiency reasons this value is 
                                                                                                  hard-coded and not currently available to be changed by the user. 

`bExp`                                  `float`     `5.0`                                         This is the second exponent in the inverse distance calc. However, for efficiency reasons this value is 
                                                                                                  hard-coded and not currently available to be changed by the user. 

`LdefFact`                              `float`     `1.0`                                         This parameter is used to determine how "far" into the field a surface deformation propagates before it is attenuated. 
                                                                                                  For small shape modifications, such as small changes to the shape of a airfoil, the default value of 1.0 is likely to be
												  sufficient. However, for much larger changes in the surface such as wing planform changes, much larger values tend to be more
												  robust. For these cases, values in the range 50-100 are common. 

`alpha`                                 `float`     `0.25`                                        This value determines how the two different exponent terms are blended. It determines how much of the higher exponent bExp
                                                                                                  term is used. Typical values are between 0.1 and 0.3. A lower value 
                                                                                                  prioritizes full blending and may result in quality reduction in the near-wall boundary layer. Higher values of alpha will 
												  tend maintain near wall quality better, but may give unacceptable skewness in the transition region between where bExp is most
												  significant to where aExp is more significant. 

`errTol`                                `float`     `0.0001`                                      This is the relative tolerance used to the fast sum approximation. A larger tolerance is faster, but may result in small 
                                                                                                  mesh imperfections away from the surface. If mesh edge lengths grow uniformly away from the body, small "errors" is the node 
												  position are not an issue. However, if the mesh has small edge lengths a great distance from the body, these imperfections may cause
												  issues and it may be required to lower the tolerance by an order of magnitude or two at the cost of more computational time. 

`evalMode`                              `str`       `fast`                                        How to compute the sums. The default which should be used at all times is 'fast'. The other option is 'exact' which is only 
                                                                                                  typically used for debugging or comparison purposes. 

`useRotations`                          `boolean`   `True`                                        Flag specifying if rotations are to be interpolated in addition to displacements. For small mesh changes it may not be necessary to
                                                                                                  interpolate rotations. However, if the surface is undergoing large changes in orientation, using rotations will help preseve 
												  boundary layer orthogonality which is generally desirable. 

`zeroCornerRotations`                   `boolean`   `True`                                        Flag specifying if rotations at sharp corners (cornerAngle defines "sharp") are to be zeroed and not contribute to the deformation. 
                                                                                                  Since the normal direction is not well defined at a corner point, including them may cause issues on some grids. 

`cornerAngle`                           `float`     `30.0`                                        The minimum deviation between surface normals surrounding a node for it to be considered a corner point. 

'bucketSize`                            `int`       `8`                                           The size of the "buckets" at the last level of the KD-tree. A large bucket size reduces the number of levels in the tree and the 
                                                                                                  overall tree size but may require more computation since a less fine granularity of leaves are available. Experiments have indicated 
												  there is little difference in run time for bucket sizes 1, 2, 4 and 8. 

======================================  ==========  ===========================================   ================================================================================================================================================================================
