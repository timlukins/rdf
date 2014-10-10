# Autogenerated with SMOP version 
# /usr/local/bin/smop decomposechange.m -o decomposechange.py
from __future__ import division
from runtime import *

def decomposechange(dS,nargout=1):
    A=matlabarray([[1,0,0,1],[1,0,0,- 1],[0,1,1,0],[0,- 1,1,0]])
    B=reshape_(dS,4,1)
    x=lscov_(A.T,B)
    x[4]=- x[4]
    E=abs_(x[1])
    S=abs_(x[2] + sqrt_(- 1) * x[3])
    R=abs_(x[4])
    return E,S,R