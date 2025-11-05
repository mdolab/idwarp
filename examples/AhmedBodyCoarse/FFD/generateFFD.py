import numpy as np


def writeFFDFile(fileName, nBlocks, nx, ny, nz, points):
    """
    Take in a set of points and write the plot 3dFile
    """

    f = open(fileName, "w")

    f.write("%d\n" % nBlocks)
    for i in range(nBlocks):
        f.write("%d %d %d " % (nx[i], ny[i], nz[i]))
    # end
    f.write("\n")
    for block in range(nBlocks):
        for k in range(nz[block]):
            for j in range(ny[block]):
                for i in range(nx[block]):
                    f.write("%f " % points[block][i, j, k, 0])
                # end
            # end
        # end
        f.write("\n")

        for k in range(nz[block]):
            for j in range(ny[block]):
                for i in range(nx[block]):
                    f.write("%f " % points[block][i, j, k, 1])
                # end
            # end
        # end
        f.write("\n")

        for k in range(nz[block]):
            for j in range(ny[block]):
                for i in range(nx[block]):
                    f.write("%f " % points[block][i, j, k, 2])
                # end
            # end
        # end
    # end
    f.close()
    return


def returnBlockPoints(corners, nx, ny, nz):
    """
    Corners needs to be 8 x 3
    """
    points = np.zeros([nx, ny, nz, 3])

    # points 1 - 4 are the iMin face
    # points 5 - 8 are the iMax face

    for idim in range(3):
        edge1 = np.linspace(corners[0][idim], corners[4][idim], nx)
        edge2 = np.linspace(corners[1][idim], corners[5][idim], nx)
        edge3 = np.linspace(corners[2][idim], corners[6][idim], nx)
        edge4 = np.linspace(corners[3][idim], corners[7][idim], nx)

        for i in range(nx):
            edge5 = np.linspace(edge1[i], edge3[i], ny)
            edge6 = np.linspace(edge2[i], edge4[i], ny)
            for j in range(ny):
                edge7 = np.linspace(edge5[j], edge6[j], nz)
                points[i, j, :, idim] = edge7
            # end
        # end
    # end

    return points


nBlocks = 4

nx = [3, 3, 3, 3]  # 4
ny = [3, 3, 3, 3]
nz = [2, 2, 2, 2]

corners = np.zeros([nBlocks, 8, 3])

corners[0, 0, :] = [-2.0, 0.0, 0.0]
corners[0, 1, :] = [-2.0, 0.0, 0.5]
corners[0, 2, :] = [-2.0, 1.0, 0.0]
corners[0, 3, :] = [-2.0, 1.0, 0.5]
corners[0, 4, :] = [0.0, 0.0, 0.0]
corners[0, 5, :] = [0.0, 0.0, 0.5]
corners[0, 6, :] = [0.0, 1.0, 0.0]
corners[0, 7, :] = [0.0, 1.0, 0.5]

corners[1, 0, :] = [0.0, 0.0, 0.0]
corners[1, 1, :] = [0.0, 0.0, 0.5]
corners[1, 2, :] = [0.0, 1.0, 0.0]
corners[1, 3, :] = [0.0, 1.0, 0.5]
corners[1, 4, :] = [0.1, 0.0, 0.0]
corners[1, 5, :] = [0.1, 0.0, 0.5]
corners[1, 6, :] = [0.1, 1.0, 0.0]
corners[1, 7, :] = [0.1, 1.0, 0.5]

corners[2, 0, :] = [0.1, 0.0, 0.0]
corners[2, 1, :] = [0.1, 0.0, 0.5]
corners[2, 2, :] = [0.1, 1.0, 0.0]
corners[2, 3, :] = [0.1, 1.0, 0.5]
corners[2, 4, :] = [1.044, 0.0, 0.0]
corners[2, 5, :] = [1.044, 0.0, 0.5]
corners[2, 6, :] = [1.044, 1.0, 0.0]
corners[2, 7, :] = [1.044, 1.0, 0.5]

corners[3, 0, :] = [1.044, 0.0, 0.0]
corners[3, 1, :] = [1.044, 0.0, 0.5]
corners[3, 2, :] = [1.044, 1.0, 0.0]
corners[3, 3, :] = [1.044, 1.0, 0.5]
corners[3, 4, :] = [5.0, 0.0, 0.0]
corners[3, 5, :] = [5.0, 0.0, 0.5]
corners[3, 6, :] = [5.0, 1.0, 0.0]
corners[3, 7, :] = [5.0, 1.0, 0.5]


# corner1 = [[-2.,0.,0.]]
# corner2 = [[5.,1.,0.5]]

# dx = (corner2[0][0]-corner1[0][0])/(nx[0]-1)
# dy = (corner2[0][1]-corner1[0][1])/(ny[0]-1)
# dz = (corner2[0][2]-corner1[0][2])/(nz[0]-1)
# print dx,dy,dz
points = []
# points.append(np.zeros([nx[0],ny[0],nz[0],3]))
for i in range(nBlocks):
    points.append(returnBlockPoints(corners[i], nx[i], ny[i], nz[i]))
# end

# for k in range(nz[0]):
#     for j in range(ny[0]):
#         for i in range(nx[0]):
#             print 'idx',i*dx
#             print 'j*dy',j*dy
#             print 'k*dz',k*dz
#             points[0][i,j,k,0] = corner1[0][0]+i*dx
#             points[0][i,j,k,1] = corner1[0][1]+j*dy
#             points[0][i,j,k,2] = corner1[0][2]+k*dz
#         # end
#     # end
# # end

# print points

fileName = "testffd.fmt"
writeFFDFile(fileName, nBlocks, nx, ny, nz, points)


nBlocks = 3

nx = [3, 5, 2]
ny = [5, 5, 5]
nz = [6, 6, 6]

corners = np.zeros([nBlocks, 8, 3])

corners[0, 0, :] = [-0.01, 0.0, 0.049]
corners[0, 1, :] = [-0.01, 0.2, 0.049]
corners[0, 2, :] = [-0.01, 0.0, 0.34]
corners[0, 3, :] = [-0.01, 0.2, 0.34]
corners[0, 4, :] = [0.1, 0.0, 0.049]
corners[0, 5, :] = [0.1, 0.2, 0.049]
corners[0, 6, :] = [0.1, 0.0, 0.34]
corners[0, 7, :] = [0.1, 0.2, 0.34]

corners[1, 0, :] = [0.1, 0.0, 0.049]
corners[1, 1, :] = [0.1, 0.2, 0.049]
corners[1, 2, :] = [0.1, 0.0, 0.34]
corners[1, 3, :] = [0.1, 0.2, 0.34]
corners[1, 4, :] = [0.862, 0.0, 0.049]
corners[1, 5, :] = [0.862, 0.2, 0.049]
corners[1, 6, :] = [0.862, 0.0, 0.34]
corners[1, 7, :] = [0.862, 0.2, 0.34]

corners[2, 0, :] = [0.862, 0.0, 0.049]
corners[2, 1, :] = [0.862, 0.2, 0.049]
corners[2, 2, :] = [0.862, 0.0, 0.34]
corners[2, 3, :] = [0.862, 0.2, 0.34]
corners[2, 4, :] = [1.05, 0.0, 0.049]
corners[2, 5, :] = [1.05, 0.2, 0.049]
corners[2, 6, :] = [1.05, 0.0, 0.21]
corners[2, 7, :] = [1.05, 0.2, 0.21]


# corner1 = [[-0.01,0.,0.049],[0.86,0.0,0.049]]
# corner2 = [[0.86,0.2,0.34],[1.05,0.2,0.21]]


points = []
for block in range(nBlocks):
    points.append(returnBlockPoints(corners[block], nx[block], ny[block], nz[block]))
    # points.append(np.zeros([nx[block],ny[block],nz[block],3]))

    # dx = (corner2[block][0]-corner1[block][0])/(nx[block]-1)
    # dy = (corner2[block][1]-corner1[block][1])/(ny[block]-1)
    # dz = (corner2[block][2]-corner1[block][2])/(nz[block]-1)
    # print dx,dy,dz

    # for k in range(nz[block]):
    #     for j in range(ny[block]):
    #         for i in range(nx[block]):
    #             print 'idx',i*dx
    #             print 'j*dy',j*dy
    #             print 'k*dz',k*dz
    #             points[block][i,j,k,0] = corner1[block][0]+i*dx
    #             points[block][i,j,k,1] = corner1[block][1]+j*dy
    #             points[block][i,j,k,2] = corner1[block][2]+k*dz
    #         # end
    #     # end
    # # end
# end

# print points

fileName = "bodyffd.fmt"
writeFFDFile(fileName, nBlocks, nx, ny, nz, points)
