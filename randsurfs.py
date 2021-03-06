# Autogenerated with SMOP version 
# /usr/local/bin/smop randsurfs.m -o randsurfs.py
from __future__ import division
from runtime import *

def randsurfs(nargout=1):
    s=0.05
    x,y=meshgrid_(arange_(- 3,3,s),arange_(- 3,3,s),nargout=2)
    z=peaks_(x,y)
    I,II=fundforms_(x,y,z,s,nargout=2)
    zi=peaks_(x,y) * 1.1
    Ii,IIi=fundforms_(x,y,zi,s,nargout=2)
    for u in arange_(1,size_(I,1)).reshape(-1):
        for v in arange_(1,size_(I,1)).reshape(-1):
            A[:,:]=Ii[u,v,:,:] - I[u,v,:,:]
            B[:,:]=IIi[u,v,:,:] - II[u,v,:,:]
            Da[u,v]=rank_(A)
            Db[u,v]=rank_(B)
    figure
    subplot_(2,1,1)
    imshow_(Da,[1,3])
    hold_('on')
    colormap_('jet')
    colorbar
    title_('A rank')
    subplot_(2,1,2)
    imshow_(Db,[1,3])
    hold_('on')
    colormap_('jet')
    colorbar
    title_('B rank')
    return
