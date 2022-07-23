import os
import numpy as np
import imageio
import math

def bilinear_interpolation(x,y,points):
		q11,q21,q22,q12,x1,y1,x2,y2 = points
		return (q11 * (x2 - x) * (y2 - y) + q21 * (x - x1) * (y2 - y) + q12 * (x2 - x) * (y - y1) + q22 * (x - x1) * (y - y1)) / ((x2 - x1) * (y2 - y1))

def s(x):
    if x >= 0:
        return 1
    else:
        return 0


def lbp(image,P=4, Radius=1):
    """LBP (Local Binary Patterns).
    Parameters
    ----------
    image --- (N, M) numpy array Graylevel image.
    P     --- Number of circularly neighbor.
    R     --- Radius of circle (spatial resolution of the operator).
    Returns
    -------
    output : (N, M) numpy array
        LBP image.
    """

    # P > 1
    # R > 0

    height = image.shape[0]
    width = image.shape[1]
    padd_image = np.zeros((height+(Radius*2),width+(Radius*2)))
    padd_image[Radius:height+1,Radius:width+1] = image
    res_image = np.zeros(image.shape)   
    
    # campute the inner matrix 
    for i in range(Radius, height):
        for j in range(Radius, width):
            lbp = 0
            for n in range(P):
                theta = 2 * math.pi * n/P
                x = i + round(Radius * math.cos((2 * math.pi * n)/P),3)
                y = j - round(-1 * Radius * math.sin((2 * math.pi * n)/P),3)
                if x.is_integer() and y.is_integer():
                    lbp = lbp + s((padd_image[int(x),int(y)]-padd_image[i,j])) * 2**n
                else:
                    # interpolation
                    points = []
                    if 0 <= theta <=  math.pi/2:
                        q11 = padd_image[i,j]
                        q21 = padd_image[i+Radius,j]
                        q22 = padd_image[i+Radius , j+Radius]
                        q12 = padd_image[i,j+Radius]
                        x1 = i 
                        y1 = j
                        x2 = i+Radius
                        y2 = j+Radius
                        points = [q11,q21,q22,q12,x1,y1,x2,y2]
                        pixel = int(bilinear_interpolation(x,y,points))
                    if math.pi/2 < theta <=  math.pi:
                        q11 = padd_image[i,j]
                        q21 = padd_image[i-Radius,j]
                        q22 = padd_image[i-Radius , j+Radius]
                        q12 = padd_image[i,j+Radius]
                        x1 = i 
                        y1 = j
                        x2 = i-Radius
                        y2 = j+Radius
                        points = [q11,q21,q22,q12,x1,y1,x2,y2]
                        pixel = int(bilinear_interpolation(x,y,points))
                    if math.pi < theta <=  3*math.pi/2:
                        q11 = padd_image[i,j]
                        q21 = padd_image[i-Radius,j]
                        q22 = padd_image[i-Radius , j-Radius]
                        q12 = padd_image[i,j-Radius]
                        x1 = i 
                        y1 = j
                        x2 = i-Radius
                        y2 = j-Radius
                        points = [q11,q21,q22,q12,x1,y1,x2,y2]
                        pixel = int(bilinear_interpolation(x,y,points))
                    if 3*math.pi/2 < theta <=  2*math.pi:
                        q11 = padd_image[i,j]
                        q21 = padd_image[i+Radius,j]
                        q22 = padd_image[i+Radius , j-Radius]
                        q12 = padd_image[i,j-Radius]
                        x1 = i 
                        y1 = j
                        x2 = i+Radius
                        y2 = j-Radius
                        points = [q11,q21,q22,q12,x1,y1,x2,y2]
                        pixel = int(bilinear_interpolation(x,y,points))
                    lbp = lbp + s((pixel-padd_image[i,j])) * 2**n 
            
            res_image[i-Radius,j-Radius] = lbp

    return res_image