import numpy as np
import matplotlib.pyplot as plt


def pixel_data(dx,dy,xsize,ysize,uvx,uvy,lam=0.21,uniform=None,t=1):
    """
    This method is make the scatter data set into pixel.
    
    """
    uv_data = np.zeros( (ysize,xsize) )
    x_ind = uvx/(dx*lam) + int(xsize/2)
    y_ind = uvy/(dy*lam) + int(ysize/2)
    if uniform == None:
        for i in range(len(y_ind)):
            uv_data[int(y_ind[i])][int(x_ind[i])] = 1 + uv_data[int(y_ind[i])][int(x_ind[i])]*np.exp( -(y_ind[i]**2+x_ind[i]**2)/t**2 )
        return uv_data
    elif uniform == "u":
        for i in range(len(y_ind)):
            uv_data[int(y_ind[i])][int(x_ind[i])] = 1 + uv_data[int(y_ind[i])][int(x_ind[i])]
        for i in range(len(uv_data)):
            for j in range(len(uv_data[0])):
                if uv_data[i][j] == 0:
                    pass
                else:
                    uv_data[i][j] = 1/uv_data[i][j]
        return uv_data





def uv_fig():
    """
    This method is to get uv figure faster.
    """
    
    return True