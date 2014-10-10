
from numpy import *
from scipy.interpolate import griddata

def quadricdata(a,b,theta,d):
    lm_sx,lm_sy=meshgrid(linspace(- 0.5,0.5,d),linspace(- 0.5,0.5,d))
    lm_sz=a*(lm_sx ** 2) + b*(lm_sy ** 2)
    ls_pts=0
    lm_nx = zeros((d*d,1))
    lm_ny = zeros((d*d,1))
    lm_nz = zeros((d*d,1))
    for x in arange(1,size(lm_sx,0)):
        for y in arange(1,size(lm_sy,1)):
            ls_pts=ls_pts + 1
            lm_nx[ls_pts]=lm_sx[x,y]
            lm_ny[ls_pts]=lm_sy[x,y]
            lm_nz[ls_pts]=lm_sz[x,y]
    R=array([[cos(theta),- sin(theta),0],[sin(theta),cos(theta),0],[0,0,1]])
    pts=array([[lm_nx],[lm_ny],[lm_nz]])
    pts=pts.T * R
    lm_nx=pts[:,1]
    lm_ny=pts[:,2]
    lm_nz=pts[:,3]
    x,y=meshgrid(linspace(- 1,1,d),linspace(- 1,1,d))
    z=griddata(array([[lm_nx],[lm_ny]]),lm_nz,(x,y))
    #z=griddata(lm_nx,lm_ny,lm_nz,x,y)
    for xi in arange(1,size(z,0)):
        for yi in arange(1,size(z,1)):
            if (isnan(z[xi,yi])):
                z[xi,yi]=0
    return x,y,z
