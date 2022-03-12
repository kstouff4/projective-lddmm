import sys
sys.path.append('../')
sys.path.append('../SurfaceTools/')

import deformSegmentations3D as ds
import numpy as np
import scipy as sp
from scipy import io
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import colors

import math
import os
import glob
import scipy.interpolate as spi
import scipy as sp
# these methods from scipy will be used for displaying some images
from scipy.linalg import eigh
from scipy.stats import norm

# Tools for Generating Paths
from numpy import linalg as LA
from itertools import combinations
from scipy.special import comb
from scipy import stats

# Tools for caching
from os import path

# Tools for saving images
import scipy.misc
from scipy import io

# Tools for Reading in File
import nibabel as nib
from nibabel import processing
from nibabel import funcs

from PIL import Image
Image.MAX_IMAGE_PIXELS
Image.MAX_IMAGE_PIXELS=1e10 # forget attack


#######################################################################
# Helper Functions
# Assuming coordinates indicate every pixel in I_orig, overlay coordVals on I_orig
def plotNew(I_orig, coordsO, coordVals, **kwargs):
    '''
    Args:
    I_orig = original image from which coordsO came
    coordsO = original pixel location coordinates in corresponding order with coordVals (assume in 2D grid)
    coordVals = value of function of interest at coordsO
    
    kwargs:
    labels = labels to label colorbar with 
    savename = filename to save image under
    
    Returns:
    I_new = function values in space of original image coordinates
    '''
    
    labels = kwargs.get('labels', None)
    savename = kwargs.get('savename',None)
    
    #I_new = np.zeros((I_orig.shape[0], I_orig.shape[1], coordVals.shape[-1])) # eliminate z coordinate
    #I_new = coordVals[:,:,0,:]
    I_new = np.squeeze(coordVals)
    I_samp = I_orig[coordsO[:,:,0],coordsO[:,:,1]] # if sample is less pixels than whole image
    
    f, ax = plt.subplots(1,2)
    ax[0].imshow(I_orig)
    ax[0].set_title("Original Image")
    im = ax[1].imshow(I_new, cmap='tab20')
    ax[1].set_title("Function Values")
    c = f.colorbar(im, ax=ax[1])
    if (labels is not None):
        c.set_ticks(list(range(len(labels))))
        c.set_ticklabels(labels)
            
    # Overlay images on top of one another
    f,ax = plt.subplots()
    ax.imshow(I_samp) # show same number of pixels as have function values for
    ax.imshow(I_new, alpha=0.5, cmap='tab20')
    ax.set_title("Sample of original image overlaid with New Function Values")
    c = f.colorbar(im, ax=ax)
    if (labels is not None):
        c.set_ticks(list(range(len(labels))))
        c.set_ticklabels(labels)
    if (savename is not None):
        f.savefig(savename)
    return I_new

def getLabelFunction(mriHeader, mriImage,brainNum,fixRes=False):
    '''
    Args:
    mriHeader = filepath of mriHeader
    mriImage = filepath of mriImage
    
    returns Lx = function of R^3 to N
    x0L, x1L, x2L = indications of each of the x at which Lx is defined (in terms of pixel locations)
    '''
    I = loadSeg(mriHeader,mriImage)
    Lx = np.asanyarray(I.dataobj)
    print("size of label functions is ")
    print(Lx.shape)
    
    aff = I.affine
    
    if (fixRes):
        d0 = np.abs(aff[0,0])
        d1 = np.abs(aff[1,1])
        d2 = np.abs(aff[2,2])
    
    else:
        # assume 1/8 mm / voxel
        d0 = 0.125
        d1 = 0.125
        d2 = 0.125

    # switch y and z
    x0L = np.arange(Lx.shape[0],dtype='float32')*d0 # scale to appropriate pixel to tissue size
    x2L = np.arange(Lx.shape[2],dtype='float32')*d2 
    x1L = np.arange(Lx.shape[1],dtype='float32')*d1 
    if (brainNum != 2 and brainNum != 1):
        x0L -= np.mean(x0L) # center at origin: 0,0 --> don't center
        x1L -= np.mean(x1L)
        x2L -= np.mean(x2L)
    if (brainNum == 7):
        print('new label function shape')
        xtemp = np.copy(x0L)
        x0L = np.copy(x2L) #np.copy(np.flip(x2L))
        x2L = xtemp
        LxNew = np.flip(Lx,2)
        Lx = np.swapaxes(LxNew,0,2)
        print(Lx.shape)
    # for brain 2 only (change the offsets to have center be at same place for 3D mri (i.e. with origin at corner):
    #x0L += (38.33999914-39.9375)
    #x1L += (29.9375-28.73999936)
    #x2L += (24.9375-23.93999946)
    X0L, X1L, X2L = np.meshgrid(x0L,x1L,x2L,indexing='ij')
    XL = np.stack((X0L,X1L,X2L),axis=-1)
    
    return Lx, XL, x0L, x1L, x2L

def loadSeg(fileHdr, fileImg):
    fileMap = nib.AnalyzeImage.make_file_map()
    fileMap['image'].fileobj = fileImg
    fileMap['header'].fileobj = fileHdr
    I = nib.AnalyzeImage.from_file_map(fileMap)
    
    # Get rid of extra dimensiosn if exist
    I = nib.funcs.squeeze_image(I)
    I.set_data_dtype(np.float64)
    print(I.affine)
    
    return I

#######################################################################
# Interpolation Functions

# Apply a function to set of points
# Note: if function does not have value assigned to desired coordinate, assign "N/A"
def applyFunction(coords, fun, d0, d1, d2=None, NN=0):
    '''
    Args:
    coords = spatial coordinates (in form of domain of fun) for which to return values for (assume 3D) (XxYxZ x 3 matrix)
    fun = function from domain of coordinates (3D) to range (any number of dimensions)
    domain = coordinates at which the function values are defined (in physical space) (XxYxZx3)

    Returns:
    coordVals = value of coords according to fun
    '''
    if (d2 is None):
        D0, D1 = np.meshgrid(d0,d1,indexing='ij')
        coordVals = interp2(D0, D1, fun, coords[:,:,0], coords[:,:,1])
    else:
        # coordVals = np.zeros((coords.shape[0], coords.shape[1], coords.shape[2], numFunValues))
        if (NN == 1):
            coordVals = interp3NN(d0,d1,d2,fun,coords[:,:,:,0], coords[:,:,:,1], coords[:,:,:,2])
        elif (NN == 2):
            coordVals = interp3Mode(d0,d1,d2,fun,coords[:,:,:,0], coords[:,:,:,1], coords[:,:,:,2])
        else:
            coordVals = interp3(d0, d1, d2, fun, coords[:,:,:,0], coords[:,:,:,1], coords[:,:,:,2])
        # Alternative interpret function that interprets based on NN rather than averaging 
    
    return coordVals

def interp2(X0, X1, I, X0s, X1s, bc='nearest'):
    ''' linear interpolation
    I want this to work for I of arbitrary dimension, interpolating across the first 2
    
    Args:
    X0 = spatial coordinates (mm) of 1st direction (used for conversion --> associated with I) (2D grid)
    X1 = spatial coordinates (mm) of 2nd direction (used for conversion --> associated with I) (2D grid)
    I = function values
    X0s = where we want to find function values
    X1s = where we want to find function values
    '''
    # convert sample points to index coords
    dx0 = X0[1,0]-X0[0,0]
    dx1 = X1[0,1]-X1[0,0]
    X0si = ((X0s - X0[0,0])/dx0).astype('float32')
    X1si = ((X1s - X1[0,0])/dx1).astype('float32')
    print('should match resolution of phi')
    print(dx0)
    print(dx1)
    
    # get fraction to next for weights
    X0si0 = np.floor(X0si).astype(int)
    X1si0 = np.floor(X1si).astype(int)
    p0 = X0si - X0si0
    p1 = X1si - X1si0
    X0si1 = X0si0+1
    X1si1 = X1si0+1
    # add necessary axes to p
    nadd = len(I.shape)-2
    for i in range(nadd):
        p0 = p0[...,None]
        p1 = p1[...,None]
    
    # boundary conditions.  This is nearest neighbor extrapolation which is usually appropriate
    if isinstance(bc,np.ndarray):
        bad00 = X0si0<0
        bad00 = np.logical_or(bad00,X0si0>I.shape[0]-1)
        bad00 = np.logical_or(bad00,X1si0<0)
        bad00 = np.logical_or(bad00,X1si0>I.shape[1]-1)
        
        bad10 = X0si1<0 # means 0 is incremented
        bad10 = np.logical_or(bad10,X0si1>I.shape[0]-1)
        bad10 = np.logical_or(bad10,X1si0<0)
        bad10 = np.logical_or(bad10,X1si0>I.shape[1]-1)
        
        bad01 = X0si0<0 # means 1 is incremented
        bad01 = np.logical_or(bad01,X0si0>I.shape[0]-1)
        bad01 = np.logical_or(bad01,X1si1<0)
        bad01 = np.logical_or(bad01,X1si1>I.shape[1]-1)
        
        bad11 = X0si1<0 # means 0 and 1 is incremented
        bad11 = np.logical_or(bad11,X0si1>I.shape[0]-1)
        bad11 = np.logical_or(bad11,X1si1<0)
        bad11 = np.logical_or(bad11,X1si1>I.shape[1]-1)
       
            
    # set boundary conditions to nearest        
    X0si0[X0si0<0] = 0
    X0si0[X0si0>I.shape[0]-1] = I.shape[0]-1
    X1si0[X1si0<0] = 0
    X1si0[X1si0>I.shape[1]-1] = I.shape[1]-1
    X0si1[X0si1<0] = 0
    X0si1[X0si1>I.shape[0]-1] = I.shape[0]-1
    X1si1[X1si1<0] = 0
    X1si1[X1si1>I.shape[1]-1] = I.shape[1]-1
    
    # vectorize (note that ravel and reshape goes across rows first)
    X0si0 = X0si0.ravel()
    X0si1 = X0si1.ravel()
    X1si0 = X1si0.ravel()
    X1si1 = X1si1.ravel()
    
    # this is if ravel would go down columns
    #X00 = X0si0 + X1si0*I.shape[0]
    #X01 = X0si0 + X1si1*I.shape[0]
    #X10 = X0si1 + X1si0*I.shape[0]
    #X11 = X0si1 + X1si1*I.shape[0]
    
    # this is if ravel goes across rows
    X00 = X0si0*I.shape[1] + X1si0
    X01 = X0si0*I.shape[1] + X1si1
    X10 = X0si1*I.shape[1] + X1si0
    X11 = X0si1*I.shape[1] + X1si1
                
    # sample four times
    # input shape
    nI = list(I.shape)
    nIravel = [nI[0]*nI[1]]
    nIravel.extend(nI[2:])
    # output shape
    n = list(X0s.shape)
    n.extend(nI[2:])
    #nravel = [n[0]*n[1]]
    #nravel.extend(n[2:])
    
    I_ = np.reshape(I,nIravel)    
    I00 = np.reshape(I_[X00],n)
    I01 = np.reshape(I_[X01],n)
    I10 = np.reshape(I_[X10],n)
    I11 = np.reshape(I_[X11],n)
    
    if isinstance(bc,np.ndarray):      
        # set out of bounds to constant
        I00[bad00] = bc
        I01[bad01] = bc
        I10[bad10] = bc
        I11[bad11] = bc
        
    # output
    return I00*((1.0-p0)*(1.0-p1)) + \
           I01*((1.0-p0)*(    p1)) + \
           I10*((    p0)*(1.0-p1)) + \
           I11*((    p0)*(    p1))

# Daniel's method of interpretation 
def interp3(X0, X1, X2, I, X0s, X1s, X2s, bc='nearest'):
    ''' linear interpolation
    I want this to work for I of arbitrary dimension, interpolating across the first 2
    
    Args:
    X0 = spatial coordinates of domain (1 x X)
    X1 = spatial coordinates of domain (1 x Y)
    X2 = spatial coordinates of domain (1 x Z)
    I = function values
    X0s = where we want to find function values (3D grid with x coord)
    X1s = where we want to find function values (3D grid with y coord)
    X2s = where we want to find function values (3D grid with z coord)
    
    returns function value (3 space)
    '''
    # convert sample points to index coords (since from matlab, should be xy only????)
    X0_1d = X0.astype('float32')
    X1_1d = X1.astype('float32')
    X2_1d = X2.astype('float32')
    
    # squeeze if first dimension is 1 and exists second 
    if len(X0.shape) > 1:
        X0_1d= np.squeeze(X0)
    if len(X1.shape) > 1:
        X1_1d = np.squeeze(X1)
    if len(X2.shape) > 1:
        X2_1d = np.squeeze(X2)
        
    # convert sample points to index coords (since from matlab, should be xy only????)
    if (X0_1d.shape[0] < 2):
        dx0 = 1 # don't know how to scale this
        X0si = (X0s - X0_1d[0])/dx0
        X0si = X0si.astype('float32')
    else:
        dx0 = X0_1d[1] - X0_1d[0]
        X0si = (X0s - X0_1d[0])/dx0
        X0si = X0si.astype('float32')
    if (X1_1d.shape[0] < 2):
        dx1 = 1 # don't know how to scale this
        X1si = (X1s - X1_1d[0])/dx1
        X1si = X1si.astype('float32')
    else:
        dx1 = X1_1d[1] - X1_1d[0]
        X1si = (X1s - X1_1d[0])/dx1
        X1si = X1si.astype('float32')
    if (X2_1d.shape[0] < 2):
        dx2 = 1 # don't know how to scale this
        X2si = (X2s - X2_1d[0])/dx2
        X2si = X2si.astype('float32')
    else:
        dx2 = X2_1d[1] - X2_1d[0]
        X2si = (X2s - X2_1d[0])/dx2
        X2si = X2si.astype('float32')
    print("Differences are: " + str(dx0) + ", " + str(dx1) + ", " + str(dx2))
    print("Should be 0.125 mm")
    
    
    # get fraction to next for weights
    X0si0 = np.floor(X0si).astype(int)
    X1si0 = np.floor(X1si).astype(int)
    X2si0 = np.floor(X2si).astype(int)
    p0 = X0si - X0si0
    p1 = X1si - X1si0
    p2 = X2si - X2si0
    X0si1 = X0si0+1
    X1si1 = X1si0+1
    X2si1 = X2si0+1
    
    # add necessary axes to p
    nadd = len(I.shape)-3
    for i in range(nadd):
        p0 = p0[...,None]
        p1 = p1[...,None]
        p2 = p2[...,None]
    
    # boundary conditions.  This is nearest neighbor extrapolation which is usually appropriate
    # 0 denotes floor being out of bounds, 1 denotes ceiling being out of bounds
    if isinstance(bc,np.ndarray):
        #floor for x, y, or z index is out of bounds
        bad000 = X0si0<0
        bad000 = np.logical_or(bad000,X0si0>I.shape[0]-1) # x0 coordinate is out of bounds
        bad000 = np.logical_or(bad000,X1si0<0) # x1 coordinate is out of bounds (small)
        bad000 = np.logical_or(bad000,X1si0>I.shape[1]-1) #x1 coordinate is out of bounds (large)
        bad000 = np.logical_or(bad000,X2si0<0)
        bad000 = np.logical_or(bad000,X2si0<I.shape[2]-1)
        
        # ceiling for x is out of bounds or floor for y,z
        bad100 = X0si1<0 # means 0 is incremented
        bad100 = np.logical_or(bad100,X0si1>I.shape[0]-1)
        bad100 = np.logical_or(bad100,X1si0<0)
        bad100 = np.logical_or(bad100,X1si0>I.shape[1]-1)
        bad100 = np.logical_or(bad100,X2si0<0)
        bad100 = np.logical_or(bad100,X2si0>I.shape[2]-1)
        
        bad010 = X0si0<0 # means 1 is incremented
        bad010 = np.logical_or(bad010,X0si0>I.shape[0]-1)
        bad010 = np.logical_or(bad010,X1si1<0)
        bad010 = np.logical_or(bad010,X1si1>I.shape[1]-1)
        bad010 = np.logical_or(bad010,X2si0<0)
        bad010 = np.logical_or(bad010,X2si0>I.shape[2]-1)
        
        bad110 = X0si1<0 # means 0 and 1 is incremented
        bad110 = np.logical_or(bad110,X0si1>I.shape[0]-1)
        bad110 = np.logical_or(bad110,X1si1<0)
        bad110 = np.logical_or(bad110,X1si1>I.shape[1]-1)
        bad110 = np.logical_or(bad110,X2si0<0)
        bad110 = np.logical_or(bad110,X2si0>I.shape[2]-1)
        
        bad001 = X0si0<0 # means 0 and 1 is incremented
        bad001 = np.logical_or(bad001,X0si0>I.shape[0]-1)
        bad001 = np.logical_or(bad001,X1si0<0)
        bad001 = np.logical_or(bad001,X1si0>I.shape[1]-1)
        bad001 = np.logical_or(bad001,X2si1<0)
        bad001 = np.logical_or(bad001,X2si1>I.shape[2]-1)
        
        bad011 = X0si0<0 # means 0 and 1 is incremented
        bad011 = np.logical_or(bad011,X0si0>I.shape[0]-1)
        bad011 = np.logical_or(bad011,X1si1<0)
        bad011 = np.logical_or(bad011,X1si1>I.shape[1]-1)
        bad011 = np.logical_or(bad011,X2si1<0)
        bad011 = np.logical_or(bad011,X2si1>I.shape[2]-1)

        bad101 = X0si1<0 # means 0 and 1 is incremented
        bad101 = np.logical_or(bad101,X0si1>I.shape[0]-1)
        bad101 = np.logical_or(bad101,X1si0<0)
        bad101 = np.logical_or(bad101,X1si0>I.shape[1]-1)
        bad101 = np.logical_or(bad101,X2si1<0)
        bad101 = np.logical_or(bad101,X2si1>I.shape[2]-1)
        
        bad111 = X0si1<0 # means 0 and 1 is incremented
        bad111 = np.logical_or(bad111,X0si1>I.shape[0]-1)
        bad111 = np.logical_or(bad111,X1si1<0)
        bad111 = np.logical_or(bad111,X1si1>I.shape[1]-1)
        bad111 = np.logical_or(bad111,X2si1<0)
        bad111 = np.logical_or(bad111,X2si1>I.shape[2]-1)
       
            
    # set boundary conditions to nearest        
    X0si0[X0si0<0] = 0
    X0si0[X0si0>I.shape[0]-1] = I.shape[0]-1
    X1si0[X1si0<0] = 0
    X1si0[X1si0>I.shape[1]-1] = I.shape[1]-1
    X2si0[X2si0<0] = 0
    X2si0[X2si0>I.shape[2]-1] = I.shape[2]-1
    X0si1[X0si1<0] = 0
    X0si1[X0si1>I.shape[0]-1] = I.shape[0]-1
    X1si1[X1si1<0] = 0
    X1si1[X1si1>I.shape[1]-1] = I.shape[1]-1
    X2si1[X2si1<0] = 0
    X2si1[X2si1>I.shape[2]-1] = I.shape[2]-1
    
    
    # vectorize (note that ravel and reshape iterate over z, then y, then x)
    # in other words, x is outer loop, y is middle, z is inner loop
    # All the floor and ceiling coordinates in a list
    X0si0 = X0si0.ravel()
    X0si1 = X0si1.ravel()
    X1si0 = X1si0.ravel()
    X1si1 = X1si1.ravel()
    X2si0 = X2si0.ravel()
    X2si1 = X2si1.ravel()
    
    # this is if ravel would go down columns
    #X00 = X0si0 + X1si0*I.shape[0]
    #X01 = X0si0 + X1si1*I.shape[0]
    #X10 = X0si1 + X1si0*I.shape[0]
    #X11 = X0si1 + X1si1*I.shape[0]
    
    # this is if ravel goes across rows
    '''
    X000 = X0si0*I.shape[1]*I.shape[2] + X1si0 + I.shape[1]*X2si0
    X010 = X0si0*I.shape[1]*I.shape[2] + X1si1 + I.shape[1]*X2si0
    X001 = X0si0*I.shape[1]*I.shape[2] + X1si0 + I.shape[1]*X2si1
    X011 = X0si0*I.shape[1]*I.shape[2] + X1si1 + I.shape[1]*X2si1
    X100 = X0si1*I.shape[1]*I.shape[2] + X1si0 + I.shape[1]*X2si0
    X110 = X0si1*I.shape[1]*I.shape[2] + X1si1 + I.shape[1]*X2si0
    X101 = X0si1*I.shape[1]*I.shape[2] + X1si0 + I.shape[1]*X2si1
    X111 = X0si1*I.shape[1]*I.shape[2] + X1si1 + I.shape[1]*X2si1
    '''
    
    dims = (I.shape[0], I.shape[1], I.shape[2])
    
    # previous part redone re: Daniel's suggestion
    '''
    X000 = X0si0*I.shape[1]*I.shape[2] + X1si0*I.shape[2] + X2si0
    X010 = X0si0*I.shape[1]*I.shape[2] + X1si1*I.shape[2] + X2si0
    X001 = X0si0*I.shape[1]*I.shape[2] + X1si0*I.shape[2] + X2si1
    X011 = X0si0*I.shape[1]*I.shape[2] + X1si1*I.shape[2] + X2si1
    X100 = X0si1*I.shape[1]*I.shape[2] + X1si0*I.shape[2] + X2si0
    X110 = X0si1*I.shape[1]*I.shape[2] + X1si1*I.shape[2] + X2si0
    X101 = X0si1*I.shape[1]*I.shape[2] + X1si0*I.shape[2] + X2si1
    X111 = X0si1*I.shape[1]*I.shape[2] + X1si1*I.shape[2] + X2si1
    '''
    
    X000 = np.ravel_multi_index([X0si0, X1si0, X2si0], dims)
    X010 = np.ravel_multi_index([X0si0, X1si1, X2si0], dims)
    X001 = np.ravel_multi_index([X0si0, X1si0, X2si1], dims)
    X011 = np.ravel_multi_index([X0si0, X1si1, X2si1], dims)
    X100 = np.ravel_multi_index([X0si1, X1si0, X2si0], dims)
    X110 = np.ravel_multi_index([X0si1, X1si1, X2si0], dims)
    X101 = np.ravel_multi_index([X0si1, X1si0, X2si1], dims)
    X111 = np.ravel_multi_index([X0si1, X1si1, X2si1], dims)
    
                
    # sample 8 times (all combinations of floor and ceiling)
    # input shape
    nI = list(I.shape)
    nIravel = [nI[0]*nI[1]*nI[2]]
    nIravel.extend(nI[3:]) # extend to accomodate all dimensions (assume that functions are at least of 1 dim)
    
    # output shape
    n = list(X0s.shape)
    n.extend(nI[3:])
    #nravel = [n[0]*n[1]]
    #nravel.extend(n[2:])
    
    I_ = np.reshape(I,nIravel).astype('float32')
    
    I000 = np.reshape(I_[X000],n).astype('float32')
    if (isinstance(bc,np.ndarray)):
        I000[bad000] = bc
    I000 = I000*((1.0-p0)*(1.0-p1)*(1.0-p2))
    
    I010 = np.reshape(I_[X010],n).astype('float32')
    if (isinstance(bc,np.ndarray)):
        I010[bad010] = bc
    I000 = I000 + ((1.0-p0)*(    p1)*(1.0-p2))*I010
    
    I010 = np.reshape(I_[X001],n).astype('float32')
    if (isinstance(bc,np.ndarray)):
        I010[bad001] = bc
    I000 = I000 + ((1.0-p0)*(1.0-p1)*(p2))*I010
    
    I010 = np.reshape(I_[X100],n).astype('float32')
    if (isinstance(bc,np.ndarray)):
        I010[bad100] = bc
    I000 = I000 + ((    p0)*(1.0-p1)*(1.0-p2))*I010
    
    I010 = np.reshape(I_[X110],n).astype('float32')
    if (isinstance(bc,np.ndarray)):
        I010[bad110] = bc
    I000 = I000 + ((    p0)*(    p1)*(1.0-p2))*I010

    I010 = np.reshape(I_[X101],n).astype('float32')
    if (isinstance(bc,np.ndarray)):
        I010[bad101] = bc
    I000 = I000 + ((    p0)*(1.0-p1)*(p2))*I010

    I010 = np.reshape(I_[X011],n).astype('float32')
    if (isinstance(bc,np.ndarray)):
        I010[bad011] = bc
    I000 = I000 + ((1.0-p0)*(    p1)*(p2))*I010

    I010 = np.reshape(I_[X111],n).astype('float32')
    if (isinstance(bc,np.ndarray)):
        I010[bad111] = bc
    I000 = I000 + ((    p0)*(    p1)*(p2))*I010

    return I000
    
    
    '''
    I000 = np.reshape(I_[X000],n).astype('float32')
    I010 = np.reshape(I_[X010],n).astype('float32')
    I001 = np.reshape(I_[X001],n).astype('float32')
    I100 = np.reshape(I_[X100],n).astype('float32')
    I110 = np.reshape(I_[X110],n).astype('float32')
    I101 = np.reshape(I_[X101],n).astype('float32')
    I011 = np.reshape(I_[X011],n).astype('float32')
    I111 = np.reshape(I_[X111],n).astype('float32')
    
    if isinstance(bc,np.ndarray):      
        # set out of bounds to constant
        I000[bad000] = bc
        I010[bad010] = bc
        I100[bad100] = bc
        I110[bad110] = bc
        I001[bad001] = bc
        I101[bad101] = bc
        I011[bad011] = bc
        I111[bad111] = bc
        
    # output (1 - p for 0s; smaller the p the closer to floor value)
    # weighted average of samples by how close index is to floor vs. ceiling 
    return I000*((1.0-p0)*(1.0-p1)*(1.0-p2)) + \
           I010*((1.0-p0)*(    p1)*(1.0-p2)) + \
           I100*((    p0)*(1.0-p1)*(1.0-p2)) + \
           I110*((    p0)*(    p1)*(1.0-p2)) + \
            I001*((1.0-p0)*(1.0-p1)*(p2)) + \
            I011*((1.0-p0)*(    p1)*(p2)) + \
            I101*((    p0)*(1.0-p1)*(p2)) + \
            I111*((    p0)*(    p1)*(p2))  
     '''

def interp3NN(X0, X1, X2, I, X0s, X1s, X2s, bc='nearest'):
    ''' NN interpolation
    I want this to work for I of arbitrary dimension, interpolating across the first 2
    
    Args:
    X0 = spatial coordinates of domain (1 x X)
    X1 = spatial coordinates of domain (1 x Y)
    X2 = spatial coordinates of domain (1 x Z)
    I = function values
    X0s = where we want to find function values (3D grid with x coord)
    X1s = where we want to find function values (3D grid with y coord)
    X2s = where we want to find function values (3D grid with z coord)
    
    returns function value (3 space)
    '''
    X0_1d = np.squeeze(X0)
    X1_1d = np.squeeze(X1)
    X2_1d = np.squeeze(X2)
    # convert sample points to index coords (since from matlab, should be xy only????)
    dx0 = X0_1d[1]- X0_1d[0]
    dx1 = X1_1d[1]- X1_1d[0]
    dx2 = X2_1d[1] - X2_1d[0]
    X0si = (X0s - X0_1d[0])/dx0
    X1si = (X1s - X1_1d[0])/dx1 
    X2si = (X2s - X2_1d[0])/dx2
    
    # debugging
    print("shapes of coords in interp3NN:")
    print(X0si.shape)
    print(X1si.shape)
    print(X2si.shape)
    
    # get fraction to next for weights
    X0si0 = np.floor(X0si).astype(int)
    X1si0 = np.floor(X1si).astype(int)
    X2si0 = np.floor(X2si).astype(int)
    p0 = X0si - X0si0
    p1 = X1si - X1si0
    p2 = X2si - X2si0
    X0si1 = X0si0+1
    X1si1 = X1si0+1
    X2si1 = X2si0+1
    
    # save lengths
    # Make sure you are within boundary
    X0si1[X0si1 >= len(np.squeeze(X0))] = len(np.squeeze(X0)) - 1
    X0si1[X0si1 < 0] = 0
    X0si0[X0si0 >= len(np.squeeze(X0))] = len(np.squeeze(X0)) - 1
    X0si0[X0si0 < 0] = 0
    X1si1[X1si1 >= len(np.squeeze(X1))] = len(np.squeeze(X1)) - 1
    X1si1[X1si1 < 0] = 0
    X1si0[X1si0 >= len(np.squeeze(X1))] = len(np.squeeze(X1)) - 1
    X1si0[X1si0 < 0] = 0
    X2si1[X2si1 >= len(np.squeeze(X2))] = len(np.squeeze(X2)) - 1
    X2si1[X2si1 < 0] = 0
    X2si0[X2si0 >= len(np.squeeze(X2))] = len(np.squeeze(X2)) - 1
    X2si0[X2si0 < 0] = 0
    
    X0si = X0si.ravel()
    X1si = X1si.ravel()
    X2si = X2si.ravel()
    X0si0 = X0si0.ravel()
    X0si1 = X0si1.ravel()
    X1si0 = X1si0.ravel()
    X1si1 = X1si1.ravel()
    X2si0 = X2si0.ravel()
    X2si1 = X2si1.ravel()
    
    X0si = (X0si <= X0si0+0.5)*X0si0 + (X0si > X0si0+0.5)*X0si1
    X1si = (X1si <= X1si0+0.5)*X1si0 + (X1si > X1si0+0.5)*X1si1
    X2si= (X2si <= X2si0+0.5)*X2si0 + (X2si > X2si0+0.5)*X2si1
    
    X0NN = X0si
    X1NN = X1si
    X2NN = X2si
    
    # vectorize (note that ravel and reshape iterate over z, then y, then x)
    # in other words, x is outer loop, y is middle, z is inner loop
    # All the floor and ceiling coordinates in a list
    '''
    X0NN = X0NN.ravel()
    X1NN = X1NN.ravel()
    X2NN = X2NN .ravel()

    p0 = p0.ravel()
    p1 = p1.ravel()
    p2 = p2.ravel()
    
    
    # Determine closest point by whether p is <= 0.5 (if less than or equal, favor lower)
    X0NN[p0 <= 0.5] = X0si0.ravel()
    X1NN[p1 <= 0.5] = X1si0.ravel()
    X2NN[p2 <= 0.5] = X2si0.ravel()
    '''

    # this is if ravel goes across rows
    XNN = X0NN*I.shape[1]*I.shape[2] + X1NN*I.shape[2] + X2NN
                
    # input shape
    nI = list(I.shape)
    nIravel = [nI[0]*nI[1]*nI[2]]
    nIravel.extend(nI[3:]) # extend to accomodate all dimensions (assume that functions are at least of 1 dim)
    
    # output shape
    n = list(X0s.shape)
    n.extend(nI[3:])
    
    I_ = np.reshape(I,nIravel)    
    INN = np.reshape(I_[XNN],n)
    
    # no need for boundary conditions??
    
    return INN

def interp3Mode(X0, X1, X2, I, X0s, X1s, X2s, bc='nearest'):
    '''
    Take a weighted mode of 8 NN
    Args:
    
    X0 = spatial coordinates of domain (1 x Y)
    X1 = spatial coordinates of domain (1 x X)
    X2 = spatial coordinates of domain (1 x Z)
    I = function values (assume certain number of unique values for labels)
    X0s = where we want to find function values (3D grid with x coord)
    X1s = where we want to find function values (3D grid with y coord)
    X2s = where we want to find function values (3D grid with z coord)
    
    returns function value (3 space)
    '''
    # convert sample points to index coords (since from matlab, should be xy only????)
    X0_1d = X0.astype('float32')
    X1_1d = X1.astype('float32')
    X2_1d = X2.astype('float32')
    
    # squeeze if first dimension is 1 and exists second 
    if len(X0.shape) > 1:
        X0_1d= np.squeeze(X0.astype('float32'))
    if len(X1.shape) > 1:
        X1_1d = np.squeeze(X1.astype('float32'))
    if len(X2.shape) > 1:
        X2_1d = np.squeeze(X2.astype('float32'))
        
    # convert sample points to index coords (since from matlab, should be xy only????)
    if (X0_1d.shape[0] < 2):
        dx0 = 1 # don't know how to scale this
        X0si = (X0s - X0_1d[0])/dx0
    else:
        dx0 = X0_1d[1] - X0_1d[0]
        X0si = (X0s - X0_1d[0])/dx0
    if (X1_1d.shape[0] < 2):
        dx1 = 1 # don't know how to scale this
        X1si = (X1s - X1_1d[0])/dx1
    else:
        dx1 = X1_1d[1] - X1_1d[0]
        X1si = (X1s - X1_1d[0])/dx1
    if (X2_1d.shape[0] < 2):
        dx2 = 1 # don't know how to scale this
        X2si = (X2s - X2_1d[0])/dx2
    else:
        dx2 = X2_1d[1] - X2_1d[0]
        X2si = (X2s - X2_1d[0])/dx2
    print("Differences are: " + str(dx0) + ", " + str(dx1) + ", " + str(dx2))
    print("Should be 0.125 mm")
    
    
    # get fraction to next for weights
    X0si0 = np.floor(X0si).astype(int)
    X1si0 = np.floor(X1si).astype(int)
    X2si0 = np.floor(X2si).astype(int)
    p0 = (X0si - X0si0).astype('float32')
    p1 = (X1si - X1si0).astype('float32')
    p2 = (X2si - X2si0).astype('float32')
    X0si1 = X0si0+1
    X1si1 = X1si0+1
    X2si1 = X2si0+1
    
    # add necessary axes to p
    nadd = len(I.shape)-3
    for i in range(nadd):
        p0 = p0[...,None]
        p1 = p1[...,None]
        p2 = p2[...,None]
    
    # boundary conditions.  This is nearest neighbor extrapolation which is usually appropriate
    # 0 denotes floor being out of bounds, 1 denotes ceiling being out of bounds
    if isinstance(bc,np.ndarray):
        #floor for x, y, or z index is out of bounds
        bad000 = X0si0<0
        bad000 = np.logical_or(bad000,X0si0>I.shape[0]-1) # x0 coordinate is out of bounds
        bad000 = np.logical_or(bad000,X1si0<0) # x1 coordinate is out of bounds (small)
        bad000 = np.logical_or(bad000,X1si0>I.shape[1]-1) #x1 coordinate is out of bounds (large)
        bad000 = np.logical_or(bad000,X2si0<0)
        bad000 = np.logical_or(bad000,X2si0<I.shape[2]-1)
        
        # ceiling for x is out of bounds or floor for y,z
        bad100 = X0si1<0 # means 0 is incremented
        bad100 = np.logical_or(bad100,X0si1>I.shape[0]-1)
        bad100 = np.logical_or(bad100,X1si0<0)
        bad100 = np.logical_or(bad100,X1si0>I.shape[1]-1)
        bad100 = np.logical_or(bad100,X2si0<0)
        bad100 = np.logical_or(bad100,X2si0>I.shape[2]-1)
        
        bad010 = X0si0<0 # means 1 is incremented
        bad010 = np.logical_or(bad010,X0si0>I.shape[0]-1)
        bad010 = np.logical_or(bad010,X1si1<0)
        bad010 = np.logical_or(bad010,X1si1>I.shape[1]-1)
        bad010 = np.logical_or(bad010,X2si0<0)
        bad010 = np.logical_or(bad010,X2si0>I.shape[2]-1)
        
        bad110 = X0si1<0 # means 0 and 1 is incremented
        bad110 = np.logical_or(bad110,X0si1>I.shape[0]-1)
        bad110 = np.logical_or(bad110,X1si1<0)
        bad110 = np.logical_or(bad110,X1si1>I.shape[1]-1)
        bad110 = np.logical_or(bad110,X2si0<0)
        bad110 = np.logical_or(bad110,X2si0>I.shape[2]-1)
        
        bad001 = X0si0<0 # means 0 and 1 is incremented
        bad001 = np.logical_or(bad001,X0si0>I.shape[0]-1)
        bad001 = np.logical_or(bad001,X1si0<0)
        bad001 = np.logical_or(bad001,X1si0>I.shape[1]-1)
        bad001 = np.logical_or(bad001,X2si1<0)
        bad001 = np.logical_or(bad001,X2si1>I.shape[2]-1)
        
        bad011 = X0si0<0 # means 0 and 1 is incremented
        bad011 = np.logical_or(bad011,X0si0>I.shape[0]-1)
        bad011 = np.logical_or(bad011,X1si1<0)
        bad011 = np.logical_or(bad011,X1si1>I.shape[1]-1)
        bad011 = np.logical_or(bad011,X2si1<0)
        bad011 = np.logical_or(bad011,X2si1>I.shape[2]-1)

        bad101 = X0si1<0 # means 0 and 1 is incremented
        bad101 = np.logical_or(bad101,X0si1>I.shape[0]-1)
        bad101 = np.logical_or(bad101,X1si0<0)
        bad101 = np.logical_or(bad101,X1si0>I.shape[1]-1)
        bad101 = np.logical_or(bad101,X2si1<0)
        bad101 = np.logical_or(bad101,X2si1>I.shape[2]-1)
        
        bad111 = X0si1<0 # means 0 and 1 is incremented
        bad111 = np.logical_or(bad111,X0si1>I.shape[0]-1)
        bad111 = np.logical_or(bad111,X1si1<0)
        bad111 = np.logical_or(bad111,X1si1>I.shape[1]-1)
        bad111 = np.logical_or(bad111,X2si1<0)
        bad111 = np.logical_or(bad111,X2si1>I.shape[2]-1)
       
            
    # set boundary conditions to nearest        
    X0si0[X0si0<0] = 0
    X0si0[X0si0>I.shape[0]-1] = I.shape[0]-1
    X1si0[X1si0<0] = 0
    X1si0[X1si0>I.shape[1]-1] = I.shape[1]-1
    X2si0[X2si0<0] = 0
    X2si0[X2si0>I.shape[2]-1] = I.shape[2]-1
    X0si1[X0si1<0] = 0
    X0si1[X0si1>I.shape[0]-1] = I.shape[0]-1
    X1si1[X1si1<0] = 0
    X1si1[X1si1>I.shape[1]-1] = I.shape[1]-1
    X2si1[X2si1<0] = 0
    X2si1[X2si1>I.shape[2]-1] = I.shape[2]-1
    
    
    # vectorize (note that ravel and reshape iterate over z, then y, then x)
    # in other words, x is outer loop, y is middle, z is inner loop
    # All the floor and ceiling coordinates in a list
    X0si0 = X0si0.ravel()
    X0si1 = X0si1.ravel()
    X1si0 = X1si0.ravel()
    X1si1 = X1si1.ravel()
    X2si0 = X2si0.ravel()
    X2si1 = X2si1.ravel()
    
    print("shapes of p's")
    print(p0.shape)
    print(p1.shape)
    print(p2.shape)
    
    p0 = p0.ravel()
    p1 = p1.ravel()
    p2 = p2.ravel()
    

    
    dims = (I.shape[0], I.shape[1], I.shape[2])
    
    X000 = np.ravel_multi_index([X0si0, X1si0, X2si0], dims)
    X010 = np.ravel_multi_index([X0si0, X1si1, X2si0], dims)
    X001 = np.ravel_multi_index([X0si0, X1si0, X2si1], dims)
    X011 = np.ravel_multi_index([X0si0, X1si1, X2si1], dims)
    X100 = np.ravel_multi_index([X0si1, X1si0, X2si0], dims)
    X110 = np.ravel_multi_index([X0si1, X1si1, X2si0], dims)
    X101 = np.ravel_multi_index([X0si1, X1si0, X2si1], dims)
    X111 = np.ravel_multi_index([X0si1, X1si1, X2si1], dims)
    
    '''
    p000 = np.ravel_multi_index([1.0-p0, 1.0-p1, 1.0-p2], dims)
    p010 = np.ravel_multi_index([1.0-p0, p1, 1.0-p2], dims)
    p001 = np.ravel_multi_index([1.0-p0, 1.0-p1, p2], dims)
    p011 = np.ravel_multi_index([1.0-p0, p1, p2], dims)
    p100 = np.ravel_multi_index([p0,1.0-p1,1.0-p2], dims)
    p110 = np.ravel_multi_index([p0,p1,1.0-p2], dims)
    p101 = np.ravel_multi_index([p0, 1.0-p1, p2], dims)
    p111 = np.ravel_multi_index([p0,p1,p2], dims)
    '''
    p000 = ((1.0-p0)*(1.0-p1)*(1.0-p2)).astype('float32')
    p010 = ((1.0-p0)*(p1)*(1.0-p2)).astype('float32')
    p001 = ((1.0-p0)*(1.0-p1)*(p2)).astype('float32')
    p011 = ((1.0-p0)*(p1)*(p2)).astype('float32')
    p100 = ((p0)*(1.0-p1)*(1.0-p2)).astype('float32')
    p110 = ((p0)*(p1)*(1.0-p2)).astype('float32')
    p101 = ((p0)*(1.0-p1)*(p2)).astype('float32')
    p111 = ((p0)*(p1)*(p2)).astype('float32')
    
                
    # sample 8 times (all combinations of floor and ceiling)
    # input shape
    nI = list(I.shape)
    nIravel = [nI[0]*nI[1]*nI[2]]
    nIravel.extend(nI[3:]) # extend to accomodate all dimensions (assume that functions are at least of 1 dim)
    
    # output shape
    n = list(X0s.shape)
    n.extend(nI[3:])
    #nravel = [n[0]*n[1]]
    #nravel.extend(n[2:])
    
    I_ = np.reshape(I,nIravel)
    numUnique = np.unique(I_)
    print("numUnique is ")
    print(numUnique)
    print("shape of p and x000")
    print(p000.shape)
    print(X000.shape)
    retLabels = np.zeros((X000.shape[0],2)) # store arg max in second and linear in first
    #print(retLabels.shape)
    #retVals = np.zeros(X000.shape) # store max linear interpolation for labels
    I_labOnly = []
    I000 = []
    I010 = []
    I001 = []
    I100 = []
    I101 = []
    I011 = []
    I111 = []
    Itot = []
    for label in numUnique:
        print("label is " + str(label))
        I_labOnly = I_ == label
        I_labOnly = I_labOnly.astype('int') # potentially store as int?
        I000 = I_labOnly[X000]
        I010 = I_labOnly[X010]
        I001 = I_labOnly[X001]
        I100 = I_labOnly[X100]
        I110 = I_labOnly[X110]
        I101 = I_labOnly[X101]
        I011 = I_labOnly[X011]
        I111 = I_labOnly[X111]
        
        '''
        Itot = I000*((1.0-p0)*(1.0-p1)*(1.0-p2)) + \
           I010*((1.0-p0)*(    p1)*(1.0-p2)) + \
           I100*((    p0)*(1.0-p1)*(1.0-p2)) + \
           I110*((    p0)*(    p1)*(1.0-p2)) + \
            I001*((1.0-p0)*(1.0-p1)*(p2)) + \
            I011*((1.0-p0)*(    p1)*(p2)) + \
            I101*((    p0)*(1.0-p1)*(p2)) + \
            I111*((    p0)*(    p1)*(p2)) 
        '''
        Itot = (I000*p000 + I010*p010 + I001*p001 + I100*p100 + I110*p110 + I101*p101 + I011*p011 + I111*p111).astype('float32')
        retLabels[Itot > retLabels[:,0],1] = label
        retLabels[:,0] = np.maximum(retLabels[:,0], Itot)
    
    retLabels = np.reshape(retLabels[:,1],n)
    
    # assume bc='nearest'

    return retLabels
    
def interp2Mode(X0, X1, I, X0s, X1s, bc='nearest'):
    '''
    Take a weighted mode of 4 NN
    Args:
    
    X0 = spatial coordinates of domain (1 x Y)
    X1 = spatial coordinates of domain (1 x X)
    X2 = spatial coordinates of domain (1 x Z)
    I = function values (assume certain number of unique values for labels)
    X0s = where we want to find function values (3D grid with x coord)
    X1s = where we want to find function values (3D grid with y coord)
    X2s = where we want to find function values (3D grid with z coord)
    
    returns function value (3 space)
    '''
    # convert sample points to index coords (since from matlab, should be xy only????)
    X0_1d = X0.astype('float32')
    X1_1d = X1.astype('float32')
    
    # squeeze if first dimension is 1 and exists second 
    if len(X0.shape) > 1:
        X0_1d= np.squeeze(X0.astype('float32'))
    if len(X1.shape) > 1:
        X1_1d = np.squeeze(X1.astype('float32'))

        
    # convert sample points to index coords (since from matlab, should be xy only????)
    if (X0_1d.shape[0] < 2):
        dx0 = 1 # don't know how to scale this
        X0si = (X0s - X0_1d[0])/dx0
    else:
        dx0 = X0_1d[1] - X0_1d[0]
        X0si = (X0s - X0_1d[0])/dx0
    if (X1_1d.shape[0] < 2):
        dx1 = 1 # don't know how to scale this
        X1si = (X1s - X1_1d[0])/dx1
    else:
        dx1 = X1_1d[1] - X1_1d[0]
        X1si = (X1s - X1_1d[0])/dx1
    print("Differences are: " + str(dx0) + ", " + str(dx1))
    print("Should be 0.125 mm")
    
    
    # get fraction to next for weights
    X0si0 = np.floor(X0si).astype(int)
    X1si0 = np.floor(X1si).astype(int)

    p0 = (X0si - X0si0).astype('float32')
    p1 = (X1si - X1si0).astype('float32')
 
    X0si1 = X0si0+1
    X1si1 = X1si0+1
    
    # add necessary axes to p
    nadd = len(I.shape)-2
    for i in range(nadd):
        p0 = p0[...,None]
        p1 = p1[...,None]

    
    # boundary conditions.  This is nearest neighbor extrapolation which is usually appropriate
    # 0 denotes floor being out of bounds, 1 denotes ceiling being out of bounds
    if isinstance(bc,np.ndarray):
        #floor for x, y, or z index is out of bounds
        bad000 = X0si0<0
        bad000 = np.logical_or(bad000,X0si0>I.shape[0]-1) # x0 coordinate is out of bounds
        bad000 = np.logical_or(bad000,X1si0<0) # x1 coordinate is out of bounds (small)
        bad000 = np.logical_or(bad000,X1si0>I.shape[1]-1) #x1 coordinate is out of bounds (large)

        
        # ceiling for x is out of bounds or floor for y,z
        bad100 = X0si1<0 # means 0 is incremented
        bad100 = np.logical_or(bad100,X0si1>I.shape[0]-1)
        bad100 = np.logical_or(bad100,X1si0<0)
        bad100 = np.logical_or(bad100,X1si0>I.shape[1]-1)
         
        bad010 = X0si0<0 # means 1 is incremented
        bad010 = np.logical_or(bad010,X0si0>I.shape[0]-1)
        bad010 = np.logical_or(bad010,X1si1<0)
        bad010 = np.logical_or(bad010,X1si1>I.shape[1]-1)
          
        bad110 = X0si1<0 # means 0 and 1 is incremented
        bad110 = np.logical_or(bad110,X0si1>I.shape[0]-1)
        bad110 = np.logical_or(bad110,X1si1<0)
        bad110 = np.logical_or(bad110,X1si1>I.shape[1]-1)
        
            
    # set boundary conditions to nearest        
    X0si0[X0si0<0] = 0
    X0si0[X0si0>I.shape[0]-1] = I.shape[0]-1
    X1si0[X1si0<0] = 0
    X1si0[X1si0>I.shape[1]-1] = I.shape[1]-1
    X0si1[X0si1<0] = 0
    X0si1[X0si1>I.shape[0]-1] = I.shape[0]-1
    X1si1[X1si1<0] = 0
    X1si1[X1si1>I.shape[1]-1] = I.shape[1]-1
    
    # vectorize (note that ravel and reshape iterate over z, then y, then x)
    # in other words, x is outer loop, y is middle, z is inner loop
    # All the floor and ceiling coordinates in a list
    X0si0 = X0si0.ravel()
    X0si1 = X0si1.ravel()
    X1si0 = X1si0.ravel()
    X1si1 = X1si1.ravel()
      
    print("shapes of p's")
    print(p0.shape)
    print(p1.shape)
    p0 = p0.ravel()
    p1 = p1.ravel()
    
    dims = (I.shape[0], I.shape[1])
    
    X000 = np.ravel_multi_index([X0si0, X1si0], dims)
    X010 = np.ravel_multi_index([X0si0, X1si1], dims)
    X100 = np.ravel_multi_index([X0si1, X1si0], dims)
    X110 = np.ravel_multi_index([X0si1, X1si1], dims)
     
    '''
    p000 = np.ravel_multi_index([1.0-p0, 1.0-p1, 1.0-p2], dims)
    p010 = np.ravel_multi_index([1.0-p0, p1, 1.0-p2], dims)
    p001 = np.ravel_multi_index([1.0-p0, 1.0-p1, p2], dims)
    p011 = np.ravel_multi_index([1.0-p0, p1, p2], dims)
    p100 = np.ravel_multi_index([p0,1.0-p1,1.0-p2], dims)
    p110 = np.ravel_multi_index([p0,p1,1.0-p2], dims)
    p101 = np.ravel_multi_index([p0, 1.0-p1, p2], dims)
    p111 = np.ravel_multi_index([p0,p1,p2], dims)
    '''
    p000 = ((1.0-p0)*(1.0-p1)).astype('float32')
    p010 = ((1.0-p0)*(p1)).astype('float32')
    p100 = ((p0)*(1.0-p1)).astype('float32')
    p110 = ((p0)*(p1)).astype('float32')
    
                
    # sample 8 times (all combinations of floor and ceiling)
    # input shape
    nI = list(I.shape)
    nIravel = [nI[0]*nI[1]]
    nIravel.extend(nI[2:]) # extend to accomodate all dimensions (assume that functions are at least of 1 dim)
    
    # output shape
    n = list(X0s.shape)
    n.extend(nI[2:])
    #nravel = [n[0]*n[1]]
    #nravel.extend(n[2:])
    
    I_ = np.reshape(I,nIravel)
    numUnique = np.unique(I_)
    print("numUnique is ")
    print(numUnique)
    print("shape of p and x000")
    print(p000.shape)
    print(X000.shape)
    retLabels = np.zeros((X000.shape[0],2)) # store arg max in second and linear in first
    #print(retLabels.shape)
    #retVals = np.zeros(X000.shape) # store max linear interpolation for labels
    I_labOnly = []
    I000 = []
    I010 = []
    I100 = []
    I110 = []
    Itot = []
    for label in numUnique:
        print("label is " + str(label))
        I_labOnly = I_ == label
        I_labOnly = I_labOnly.astype('int') # potentially store as int?
        I000 = I_labOnly[X000]
        I010 = I_labOnly[X010]
        I100 = I_labOnly[X100]
        I110 = I_labOnly[X110]
        
        '''
        Itot = I000*((1.0-p0)*(1.0-p1)*(1.0-p2)) + \
           I010*((1.0-p0)*(    p1)*(1.0-p2)) + \
           I100*((    p0)*(1.0-p1)*(1.0-p2)) + \
           I110*((    p0)*(    p1)*(1.0-p2)) + \
            I001*((1.0-p0)*(1.0-p1)*(p2)) + \
            I011*((1.0-p0)*(    p1)*(p2)) + \
            I101*((    p0)*(1.0-p1)*(p2)) + \
            I111*((    p0)*(    p1)*(p2)) 
        '''
        Itot = (I000*p000 + I010*p010 + I100*p100 + I110*p110).astype('float32')
        retLabels[Itot > retLabels[:,0],1] = label
        retLabels[:,0] = np.maximum(retLabels[:,0], Itot)
    
    retLabels = np.reshape(retLabels[:,1],n)
    
    # assume bc='nearest'

    return retLabels
    

#######################################################################
# Slice info 
def getInfo(brainNum, stainNum, blockNum, sliceNum, locationNum,suffix='_newmai0803.npz',redone='redone_11222',suffixSave=''):
    '''
    Args: 
    sliceNum = 0 indexed out of slices within same stain
    locationNum = as saved in file
    points = list of 3D coordinates to interrogate 
    manual if to add to output name (i.e. if alignment is based on manual adjustment rather than pipeline)
    
    '''
    stainString = 'Tau'
    stain = 'PHF-1'
    if (stainNum == 0):
        stainString = 'Amyloid'
        stain = '6E10'
    # Brain and block info (coordinates will have been saved together)
    labels = []
    if (brainNum == 2):
        labels = ['background','HATA','Alveus','Amygdala','CA1','CA2','CA3','ERC-extension','ERC','DG-granular','DG-hilus','DG-molecular','Parasubiculum','Presubiculum','Subiculum'] # 1/12/22
    
    elif (brainNum == 5):
        labels = ['background','alveus','ca1','ca2','ca3','endfolial','erc','granular','HATA','hilus','molecular','parasubiculum','presubiculum','subiculum','claire_extension','amygdala']
    
    elif (brainNum == 3):
        labels = ['background','alveus','ca2','ca3','endfolial','granular','hilus','molecular','ca1','subiculum','erc','HATA','presubiculum','parasubiculum','amygdala']
    # Individual Slice info 
    #sliceNum = 5 # out of total, 0 indexed 
    #locationNum = 10 # as saved in original histology name 
    saveName = '/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/2DHisto_to_MRI_coords_block' + str(blockNum+1) + '_stain' + str(stainNum) + '_slice' + str(sliceNum) + suffix
    
    params = np.load(saveName)
    coordsT_comb = params['coordsT_comb'] # check on this 
    coordsH = params['coordsH']
    coords_2DMRI = params['coords_2DMRI']
    coordsT = params['coordsT']
    coordsT_mai = params['coordsT_mai']
    res0 = params['res0']
    res1 = params['res1']
    print('res0 is ' + str(res0))
    coordsH_pix = params['coordsH_pix']
    params.close()

    ADHip = '/AD_Hip' + str(blockNum+1)
    if (brainNum == 2):
        if blockNum == 0:
            # plot sequence of grids
            Ihisto = '/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/histology/down_000/AD_Hip' + str(blockNum+1) + '/Tau/BRC2614_AD_Block' + str(blockNum+1) + '_PHF-1_Location_' + str(locationNum) + '_corrected.tif'
        else:
            Ihisto = '/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/histology/down_000/AD_Hip' + str(blockNum+1) + '/AD_Hip' + str(blockNum+1) + '/Tau/BRC2614_AD_Block' + str(blockNum+1) + '_PHF-1_Location_' + str(locationNum) + '_corrected.tif'

        if (stainNum == 0):
            Ihisto = '/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain' + str(brainNum) + '/histology/down_000/AD_Hip' + str(blockNum+1)+ADHip + '/Amyloid/BRC2614 AD Block' + str(blockNum+1) + ' 6E10 Location ' + str(locationNum) + '_corrected.tif'
    elif (brainNum == 5 or brainNum == 3):
        Ihisto = '/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain' + str(brainNum) + '/histology/down_000' + ADHip + '/' + stainString + '/Brain ' + str(brainNum) + '-Block ' + str(blockNum+1) + ' L' + str(locationNum) + ' ' + stain + '_crop_filterNS.tif'
        
    # Get and apply Eileen's function
    mriHeader = '/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain' + str(brainNum) + '/3DSegmentations/' + redone + '/allMerge_smoothe.hdr'
    mriImage = '/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain' + str(brainNum) + '/3DSegmentations/' + redone + '/allMerge_smoothe.img'

    Lx, XL, x0L, x1L, x2L = getLabelFunction(mriHeader, mriImage,brainNum)
    L_histo_Mode = applyFunction(coordsT_comb,Lx,x0L,x1L,x2L,2) # Mode

    saveNameLab = '/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/2DHisto_to_MRI_coords_block' + str(blockNum+1) + '_stain' + str(stainNum) + '_loc' + str(locationNum) + '_labelIm_' + redone + suffixSave + '.png'
    
    L_histo_image_Mode = plotNew(plt.imread(Ihisto), coordsH_pix.astype(int), L_histo_Mode, labels=labels,savename=saveNameLab)
    
    # save image labels as npz rather than image
    saveName = '/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/2DHisto_to_MRI_coords_block' + str(blockNum+1) + '_stain' + str(stainNum) + '_loc' + str(locationNum) + '_Newlabels_' + redone + suffixSave + '.npz'
    np.savez(saveName,L_histo_Mode=L_histo_Mode) #check these are correct things to save 
   
    return
##################################################################
# Histo to MRI Functions

# Transformation From File Functions (i.e. read in mat or .py)

# get transformations and relevant parameters from Daniel's alignment program
# returns all slice transformations, use selectSliceTrans() to choose correct one
# first transformation is for the deformation vs. second one contains the XM, YM, ZM coordinates
def get3DTrans(matFileParams,matFileIM):
    params = sp.io.loadmat(matFileParams, appendmat=False, struct_as_record=False)
    
    ASi = np.asarray(params['Bsave']) # 4x4 matrix for each stain and slice index

    phiiX = np.asarray(params['phisaveX'])
    phiiY = np.asarray(params['phisaveY'])
    phiiZ = np.asarray(params['phisaveZ'])
    
    coords = sp.io.loadmat(matFileIM, appendmat=False, struct_as_record=False)
    xM = np.asarray(coords['xM']).astype('float32')
    yM = np.asarray(coords['yM']).astype('float32')
    zM = np.asarray(coords['zM']).astype('float32')
    print("shape of coords:")
    print(xM.shape)
    print(yM.shape)
    print(zM.shape)
    

    
    return ASi,xM,yM,zM,phiiX, phiiY, phiiZ
    #return Asi,X0D,X1D,X2D,phiiX,phiiY,phiiZ

# TO DO: might want to switch rows and columns 
# xM, yM, zM = all the same for every slice because in 3D MRI space 
def selectSliceTrans(ASi,phiiX, phiiY, phiiZ,stainI,sliceI):
    ASi_si = ASi[stainI,sliceI]

    phiiX_si = phiiX[stainI,sliceI].astype('float32')
    phiiY_si = phiiY[stainI,sliceI].astype('float32')
    phiiZ_si = phiiZ[stainI,sliceI].astype('float32')
    '''
    phiiX_si = phiiX[stainI,sliceI]
    phiiY_si = phiiY[stainI,sliceI]
    phiiZ_si = phiiZ[stainI,sliceI]
    phii_si = np.stack((phiiX_si,phiiY_si,phiiZ_si),axis=-1)
    '''
    # Assume phiis are the same for each stain and slice
    return ASi_si,phiiX,phiiY,phiiZ

# Assume 3D or 4D matrix gotten from matlab
# Put into equivalent python coordinates by reversing x and y (A,B), (A',B'), (A'', B'')
# (A_,B_), (A__,B__), (A___,B___)
# put X,Y,andZ coordinates together 
def transposeXandY(A, x0, x1, x2, phiiX, phiiY, phiiZ):
    # A = 4x4 double --> transpose 2x2 submatrix
    print("transpose x and y shape coords")
    print(phiiX.shape)
    print(x0.shape)
    
    Aret = np.copy(A)
    Aret[0,0] = A[1,1]
    Aret[1,1] = A[0,0]
    Aret[0,1] = A[1,0]
    Aret[1,0] = A[0,1]
    Aret[0,2] = A[1,2]
    Aret[1,2] = A[0,2]
    Aret[0,3] = A[1,3]
    Aret[1,3] = A[0,3]
    Aret[2,0] = A[2,1]
    Aret[2,1] = A[2,0]
    
    # don't change the x,y,z
    x1ret = x0
    x0ret = x1
    
    if (x2 is None):
        phiRet = np.zeros((x1.shape[1], x0.shape[1], 1, 3))
        phiRet[:,:,:,0] = phiiY[:,:,None].astype('float32')
        phiRet[:,:,:,1] = phiiX[:,:,None].astype('float32')
        phiRet[:,:,:,2] = phiiZ[:,:,None].astype('float32')
    else:
        phiRet = np.zeros((x1.shape[1],x0.shape[1],x2.shape[1],3)) # no need to transpose
        phiRet[:,:,:,0] = phiiY.astype('float32')
        phiRet[:,:,:,1] = phiiX.astype('float32')
        phiRet[:,:,:,2] = phiiZ.astype('float32')

    return Aret, x0ret, x1ret, x2, phiRet

# get transformations and relevant parameters from Katie's alignment program
# only 1 slice transformed here 
def get2DTrans(pyFile):
    params = np.load(pyFile)
    X0K = params['X0I'] # check on this 
    X1K = params['X1I']
    psi = params['phi']
    A2 = params['A']
    params.close()
    
    return X0K,X1K,psi,A2

def comBlocks(matFile, blockNum, coords, Ashift):
    params = sp.io.loadmat(matFile, appendmat=False, struct_as_record=False)
    coordsRet = np.copy(coords)
    
    A = np.asarray(params['A']) # array of 4x4 matrices for each block
    # Transpose A matrix
    
    blocks = A.shape[1] # should be 1 x 3 x 4 x 4
    print("affine matrix as is")
    Ablock = A[0,blockNum]
    print(Ablock)
    
    Aret = np.copy(Ablock)
    Aret[0,0] = Ablock[1,1]
    Aret[1,1] = Ablock[0,0]
    Aret[0,1] = Ablock[1,0]
    Aret[1,0] = Ablock[0,1]
    Aret[0,2] = Ablock[1,2]
    Aret[1,2] = Ablock[0,2]
    Aret[0,3] = Ablock[1,3]
    Aret[1,3] = Ablock[0,3]
    Aret[2,0] = Ablock[2,1]
    Aret[2,1] = Ablock[2,0]
    Ablock = Aret
    print("affine transposed")
    print(Ablock)
    if (Ashift is not None):
        Ablock = Ablock@Ashift
        
    print("Affine matrix used")
    print(Ablock)
    coordsRet[:,:,:,0] = Ablock[0,0]*coords[:,:,:,0] + Ablock[0,1]*coords[:,:,:,1] + Ablock[0,2]*coords[:,:,:,2] + Ablock[0,3]
    coordsRet[:,:,:,1] = Ablock[1,0]*coords[:,:,:,0] + Ablock[1,1]*coords[:,:,:,1] + Ablock[1,2]*coords[:,:,:,2] + Ablock[1,3]
    coordsRet[:,:,:,2] = Ablock[2,0]*coords[:,:,:,0] + Ablock[2,1]*coords[:,:,:,1] + Ablock[2,2]*coords[:,:,:,2] + Ablock[2,3]
           
    return coordsRet

def transToMai(matFile, coords, inverse=-1):
    '''
    Args: inverse: if = 1, assume coords in mai space and translate to 3D MRI space
    
    '''
    params = sp.io.loadmat(matFile, appendmat=False, struct_as_record=False)
    coordsRet = np.copy(coords)
    coordsMicron = np.copy(coords)
    #coordsMicron = 1000.0*coordsMicron
    
    A = np.asarray(params['A']) # 4x4 array only of the transformation to apply (takes microns to microns)
    
    Aret = np.copy(A)
    Aret[0,0] = A[1,1]
    Aret[1,1] = A[0,0]
    Aret[0,1] = A[1,0]
    Aret[1,0] = A[0,1]
    Aret[0,2] = A[1,2]
    Aret[1,2] = A[0,2]
    Aret[0,3] = A[1,3]
    Aret[1,3] = A[0,3]
    Aret[2,0] = A[2,1]
    Aret[2,1] = A[2,0]
    A = Aret
    if (inverse > 0):
        A = np.linalg.inv(Aret)
    #A = A/1000.0 # transform to taking and returning mm as given in Mai Atlas (want output to be in mm, not microns)
    coordsRet[:,:,:,0] = A[0,0]*coordsMicron[:,:,:,0] + A[0,1]*coordsMicron[:,:,:,1] + A[0,2]*coordsMicron[:,:,:,2] + A[0,3]
    coordsRet[:,:,:,1] = A[1,0]*coordsMicron[:,:,:,0] + A[1,1]*coordsMicron[:,:,:,1] + A[1,2]*coordsMicron[:,:,:,2] + A[1,3]
    coordsRet[:,:,:,2] = A[2,0]*coordsMicron[:,:,:,0] + A[2,1]*coordsMicron[:,:,:,1] + A[2,2]*coordsMicron[:,:,:,2] + A[2,3]
    
    return coordsRet # return in mm

# Transform set of points in 2D histology to 3D MRI
# All transformations given in form of matrices
# single letters assumed to be the same transformation applied at each position (rigid) vs. greek letters denote nonrigid
def histoToMRIPoints(coordsH, psi, A2, ASi, phii, z, res0, res1, **kwargs):
    '''
    Args:
    coordsH = histology coordinates to sample (i.e. transform to 3D MRI space) (in space, grid format: X x Y x 2)

    Katie's functions
    psi = transformation from histology coordinates to intermediate (nonrigid deform)
    A2 = affine transformation from intermediate to 2D mri
   
    Daniel's functions
    ASi = inverse of 2D rigid(translation), scaling, and affine (assumes HP units to mm and takes in 3D HP coordinate)
    A1i = inverse of 2D rigid (translation) from 2D mri to 2D pre-mri
    Si = inverse of scaling and slicing from 2D pre-mri to 3D slice position
    Ai = inverse of 3D affine from slice position to intermediate 3D mri
    phii = inverse of nonrigid deformation from intermediate to 3D mri blocks
    Bi = inverse of affine transformtion from individual blocks to combined
    z = slice number to generate corresponding z coordinate
    
    res0 = resolution of histology in 0 direction (i.e. mm/HP)
    res1 = resolution of histology in 1 direction (i.e. mm/HP)
    
    kwargs
    x0K, x1K = mm at which Katie's functions are defined in domain (used for interpolation); if not provided, assume coordsH is equivalent
    x0D, x1D = mm at which Daniel's functions are defined in domain (used for interpolation); assume in the form of 3D
    X0M,X1M = sampling coordinates of 2D MRI to transform
    Returns:
    coordsM = corresponding 3D mri coordinates
    weights = for incorporating in varifold model (should still be multiplied by 0.00208 mm^2)
    '''
    X0m = kwargs.get('X0m',None)
    X1m = kwargs.get('X1m',None)

    x0D = kwargs.get('x0D', None)
    x1D = kwargs.get('x1D', None) # might need more than one daniel domains / conversions 
    x2D = kwargs.get('x2D',None)


    if (X0m is None or X1m is None):
        #coordsM = np.zeros((coordsH.shape[0],coordsH.shape[1],1,3))
        #W2 = np.zeros((coordsH.shape[0],coordsH.shape[1])) # 2D jacobian weights
        #W3 = np.zeros((coordsH.shape[0],coordsH.shape[1])) # 3D jacobian weights
        XH = np.copy(coordsH)
        '''
        # Turn coordinates into grid (assume coordsH are given in mm)
        x0H = coordsH[:,0]
        x1H = coordsH[:,1]
        X0H,X1H = np.meshgrid(x0H,x1H,indexing='ij') # starting histology coordinates of physical space 
        XH = np.stack((X0H,X1H),axis=-1)
        '''
        # Deformation conversions 
        X0K = kwargs.get('X0K',XH[:,:,0])
        X1K = kwargs.get('X1K', XH[:,:,1])
        # calculate sampling intervals (unnecessary?)
        '''
        dx0H = x0H[1] - x0H[0]
        dx1H = x1H[1] - x1H[0]
        '''

        # Apply psi (mm/pix)
        #X0K,X1K = np.meshgrid(x0K,x1K,indexing='ij') # starting I coordinates of physical space 
        XK = np.stack((X0K,X1K),axis=-1)
        X0s = interp2(X0K, X1K, psi[...,0]-X0K, XH[:,:,0], XH[:,:,1]) + XH[...,0] #returns X x Y x 2 displacement
        X1s = interp2(X0K,X1K, psi[...,1]-X1K, XH[:,:,0], XH[:,:,1]) + XH[...,1]
        #w0,w1 = np.gradient(Xs,res0,res1,axis=(0,1)) # jacobian factor (here, assumes histology is evenly spaced at given resolution in both directions)
        #W2 = np.abs(w0[...,0]*w1[...,1] - w0[...,1]*w1[...,0]) # jacobian determinant

        #X0s = Xs[:,:,0]
        #X1s = Xs[:,:,1]
        Xs = np.stack((X0s,X1s),axis=-1)

        # Apply A2 (2 x 3 matrix) (mm/pix)
        X0m = A2[0,0]*X0s + A2[0,1]*X1s + A2[0,2]
        X1m = A2[1,0]*X0s + A2[1,1]*X1s + A2[1,2]
        #W2 = W2*np.abs(np.linalg.det(A2))# multiply by determinant of affine matrix 
    
    X_2DMRI = np.stack((X0m,X1m),axis=-1)
    
    # add z coordinate
    X0m = X0m[:,:,None]
    X1m = X1m[:,:,None]
    X2m = np.zeros(X0m.shape)
    X2m[:,:,0] = z # assume z is equal to slice number - mean(1:numSlices)
    print("z is " + str(z))
    print("X0m shape")
    print(X0m.shape)
    print(X0m[0,0])
    print(X0m[-1,-1])
    print(X1m[0,0])
    print(X1m[-1,-1])
    
    # Apply conversion factor from mm to HP or sliceI per pixel 
    #X0m = X0m*(1/0.00208)
    #X1m = X1m*(1/0.00208)
    X0mc = X0m*(1.0/res0) # turn into histology pixels
    X1mc = X1m*(1.0/res1) # turn into histology pixels 
    X2mc = X2m # already in slice number 
    
    #### KATIE NEW FIX 3/2/22 ####
    '''
        Matlab code assumes reversed rows and columns for 2D coordinates??
    '''
    X0m = X0mc
    X2m = X2mc
    X1m = X1mc
    
    print("ASi is ")
    print(ASi)
    
    # Apply combined affine transformations (B matrix) --> returns to mm/pixel
    X0m_n = ASi[0,0]*X0m + ASi[0,1]*X1m + ASi[0,2]*X2m + ASi[0,3]
    X1m_n = ASi[1 ,0]*X0m + ASi[1,1]*X1m + ASi[1,2]*X2m + ASi[1,3]
    X2m_n = ASi[2,0]*X0m + ASi[2,1]*X1m + ASi[2,2]*X2m + ASi[2,3]
    X0m = X0m_n
    X1m = X1m_n
    X2m = X2m_n
    Xm = np.stack((X0m,X1m,X2m),axis=-1)
    print("coordinates after 3D affine deformation")
    print(Xm[0,0,0,:])
    print(Xm[-1,-1,-1,:])
       
    # Apply phii
    if (x0D is None or x1D is None or x2D is None):
        print("Cannot carry out 3D deformation. Please provide domain.")
        Xm = np.zeros((X0m.shape[0],X0m.shape[1],X0m.shape[2],3))
    else:
        X0D, X1D, X2D = np.meshgrid(x0D,x1D,x2D,indexing='ij')
        XD = np.stack((X0D,X1D,X2D),axis=-1)
        print("phii based on")
        print(XD[0,0,0,:])
        print(XD[-1,-1,-1,:])
        Xm = interp3(x0D, x1D, x2D, phii-XD, X0m, X1m, X2m) + Xm # Assumes regular sampling generated the phii?
        #Xm0 = interp3(x0D,x1D,x2D,phii[...,0]-X0D,X0m,X1m,X2m) + X0m
        #Xm1 = interp3(x0D,x1D,x2D,phii[...,1]-X1D,X0m,X1m,X2m) + X1m
        #Xm2 = interp3(x0D,x1D,x2D,phii[...,2]-X2D,X0m,X1m,X2m) + X2m
        #Xm = np.stack((Xm0,Xm1,Xm2),axis=-1)
        
    #Xmrav = np.reshape(Xm,(np.prod(Xm.shape[:-1]),Xm.shape[-1]))
    #w0 = np.gradient(Xmrav,np.ravel(X0m)) # might need to add eps if Xm coordinates not different
    #w1 = np.gradient(Xmrav,np.ravel(X1m))
    #w2 = np.gradient(Xmrav,np.ravel(X2m))
    #w1[...,1]*w2[...,2] -  #finish 3d gradient
    
    # Save all of the intermediate coordinates for testing
    print("coordinates after 3D phi deformation")
    print(Xm[0,0,0,:])
    print(Xm[-1,-1,-1,:])
    print("shape of coords are:")
    print(X_2DMRI.shape)
    
    return X_2DMRI, Xm

# Generate list of Coordinates from Image as samples at every pixel (assume 2D image)
# meshed grid (assume interest in sampling at regular interval)

def generateCoords2(I, res0, res1):
    '''
    Args:
    imagename = imagename that wish to get information on
    res = mm/pixel in given image
    
    Returns:
    coordsP = pixel locations arranged in order of 2D array
    coordsS = spatial locations arranged in order of 2D array (corresponding to pixel loc ordering)
    '''
    
    # Old Version
    '''
    numLocs = I.shape[0]*I.shape[1]
    coordsP = np.array([(x,y) for x in range(I.shape[0]) for y in range(I.shape[1])])
    coordsS = np.copy(coordsP)
    coordsS[:,0] = coordsS[:,0]*res0
    coordsS[:,1] = coordsS[:,1]*res1
    '''
    I = plt.imread(I)
    # New Version (don't shift the pixel locations)
    x0 = np.arange(I.shape[0], dtype='float32')
    x0s = x0*res0
    x0s = x0s - np.mean(x0s)
    x1 = np.arange(I.shape[1], dtype='float32')
    x1s = x1*res1
    x1s = x1s - np.mean(x1s)

    
    X0,X1 = np.meshgrid(x0,x1,indexing='ij')
    coordsP = np.stack((X0,X1),axis=-1).astype('float32')
    X0s,X1s = np.meshgrid(x0s,x1s,indexing='ij')
    coordsS= np.stack((X0s,X1s),axis=-1).astype('float32') 
    
    return coordsP, coordsS

##################################################################
# Wrapper Functions

# Run code for block and save in folder
# This works off of the premise that there only exist MRI slices per slice and not per stain 
def runPipelineForBlock(brainNum, blockNum, stainNum, meanSliceNum, Ihisto, mriParams, histTrans, mriTrans, mriCoords,X0m=None,X1m=None,saveName=None,lowBound=-1,highBound=20,numSlices=None):
    '''
    Args:
    brainNum = 1 index brain Num
    blockNum = 0 index brain Num
    stainNum = 0 index stain (i.e. Amyloid, LFB, Tau for brain 2)
    meanSliceNum = meanSliceNum used to shift coordinates in Z direction by
    Ihisto = filenames to original resolution of histology (one per slice)
    mriParams = parameters for whole block 
    histTrans = filenames (one per slice) for histology transformation parameters
    mriTrans = transformations for whole block 
    mriCoords = coordinates that diffeomorphism in mriTrans is defined at
    '''

    XMstring = ''
    if (X0m is not None or X1m is not None):
        XMstring = 'XM'
    params = sp.io.loadmat(mriParams, appendmat=False, struct_as_record=False)    
    sxy = np.asarray(params['sxy'])
    sxy = sxy[0]

    res0 = 1.0/sxy # histology resolution
    res1 = 1.0/sxy # histology resolution
    
    ASi,xM,yM,zM,phiiX,phiiY,phiiZ = get3DTrans(mriTrans,mriCoords)
    print(ASi.shape)
    print('xM shape')
    print(xM.shape)
    print('phiiX shape, from MRI trans')
    print(phiiX.shape)
    
    if (numSlices is None):
        numSlices = len(Ihisto)
        
    for sl in range(numSlices):
        # Buffer for only doing certain slices 
        
        if (sl <= lowBound or sl >= highBound):
            print("skipping slice " + str(sl))
            continue

        coordsH_pix,coordsH_space = generateCoords2(Ihisto[sl],res0,res1)
        print('Histology coords')
        print(coordsH_pix.shape)
        z = (sl + 1) - meanSliceNum # 1 indexed 

        ASi_s,phiiX_s, phiiY_s, phiiZ_s = selectSliceTrans(ASi,phiiX, phiiY, phiiZ,0,0) # assume there is only one set of slices and transformations
        ASi_s, x0M, x1M, x2M, phii = transposeXandY(ASi_s,xM,yM,zM,phiiX_s, phiiY_s, phiiZ_s)

        # Carryout transformation of space (within 1 block)
        if (X0m is None or X1m is None):
            X0K, X1K, psi, A2 = get2DTrans(histTrans[sl])
            kwargs = {"X0K": X0K, "X1K": X1K, "x0D": x0M, "x1D": x1M, "x2D": x2M}
            coords_2DMRI, coordsT = histoToMRIPoints(coordsH_space, psi, A2, ASi_s, phii, z, res0, res1, **kwargs)
        else:
            psi = np.stack((X0m,X1m),axis=-1)
            X0K = X0m
            X1K = X1m
            A2 = np.eye(3)
            kwargs = {"X0K": X0K, "X1K": X1K, "x0D": x0M, "x1D": x1M, "x2D": x2M, "X0m":X0m, "X1m":X1m}
            coords_2DMRI, coordsT = histoToMRIPoints(coordsH_space, psi, A2, ASi_s, phii, z, res0, res1, **kwargs)
        print("shape of phi coords")
        print(X0K.shape)
        print("shape of desired coords")
        print(coordsH_space.shape)

        # Transform coordinates into space of combined block MRI
        # Brain1 has multiple files currently (7/3/20); amyHippo = block1, hippo = block2, caudate=block3
        if (brainNum == 2):
            blockCombMat = '/cis/home/dtward/Documents/exvivohuman_11T/more_blocks/Brain2/A_may10_touch_up_block_1.mat' # filepath for .mat saving block combinations
        elif (brainNum == 3 or brainNum == 4 or brainNum == 5 or brainNum == 7):
            blockCombMat = '/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/Params/singleToCombMRI.mat'
        else:
            print("Cannot combine blocks")
            blockCombMat = ''
        if (brainNum == 2):
            blockString = str(blockNum + 1)
            Ashift_mat = np.load('/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain2/block' + blockString + '_Ashift.npz')
            Ashift = Ashift_mat['Ashift']
            coordsT_comb = comBlocks(blockCombMat, blockNum, coordsT, Ashift)
        else:
            coordsT_comb = comBlocks(blockCombMat,blockNum,coordsT,np.eye(4)) # need to generate shift matrix for other brains

        # Transform coordinates into MAI atlas space with additional rigid transformation
        if (brainNum == 2):
            maiAtlasMat = '/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/Params/A_toMai_0803.mat'
        elif (brainNum == 5):
             maiAtlasMat = '/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/Params/A_toMai_1127.mat'
        else:
            print("cannot combine mai atlas coordinates yet")
            maiAtlasMat = np.eye(4)
        coordsT_mai = transToMai(maiAtlasMat,coordsT_comb)

        # generate 3D original coordinates 
        coordsH_space0 = coordsH_space[:,:,0].astype('float32')
        coordsH_space1 = coordsH_space[:,:,1].astype('float32')
        coordsH_space2 = np.zeros(coordsH_space0.shape)
        coordsH_space2[:,:] = z
        coordsH_space0 = coordsH_space0[...,None].astype('float32')
        coordsH_space1 = coordsH_space1[...,None].astype('float32')
        coordsH_space2 = coordsH_space2[...,None].astype('float32')
        
        coordsH = np.stack((coordsH_space0,coordsH_space1,coordsH_space2), axis=-1)
        
        # Save Coordinates
        if (saveName is None):
            saveNameNew = '/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/2DHisto_to_MRI_coords_block' + str(blockNum+1) + '_stain' + str(stainNum) + '_slice' + str(sl) + XMstring + '.npz'
        else:
            saveNameNew = '/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/2DHisto_to_MRI_coords_block' + str(blockNum+1) + '_stain' + str(stainNum) + '_slice' + str(sl) + saveName + XMstring + '.npz'
        np.savez(saveNameNew,coordsH=coordsH,coordsT=coordsT,coordsT_comb=coordsT_comb,coordsT_mai=coordsT_mai,coordsH_pix=coordsH_pix, coords_2DMRI=coords_2DMRI, res0=res0, res1=res1) #check these are correct things to save 
    
    return

def getComposition(brainNum,blockNum,stain,dateNew,mriSuffix):
    '''
    Use with script runCompose.sh
    brainNum = 1 based digit
    blockNum = 1 based digit
    stain = 'Tau' or 'Amyloid'
    mriSuffix = '' or '_pil' (before .tif)
    dateNew = folder under output figs where 2D-2D transformations are stored (e.g. '020222')
    '''
    
    # Step 1 
    # Run full pipeline on all slices of 1 block (assume same stain); (brainnum-1,blocknum-1,stain)
    sliceNums = np.zeros((7,4,3)).astype(int)
    sliceNums[1,0,2] = 15
    sliceNums[1,0,0] = 15
    sliceNums[1,1,2] = 15
    sliceNums[1,1,0] = 15
    sliceNums[1,2,2] = 7
    sliceNums[1,2,0] = 7
    sliceNums[2,0,2] = 13
    sliceNums[2,1,2] = 12
    sliceNums[2,2,2] = 13
    sliceNums[2,3,2] = 11
    sliceNums[4,0,2] = 13
    sliceNums[4,1,2] = 14
    sliceNums[4,2,2] = 8

    stainNum = 2 # 0 indexed: Amyloid, KFB, Tau
    stainString = "Tau"
    if (stain == 'Amyloid'):
        stainString = "Amyloid"
        stainNum = 0

    stain = ['6E10','LFB H&E','PHF-1']
    sliceName = ['1','5','10','15','20','25','30','35','40','45','50','55','60','65','70','75']
    
    if (brainNum == 5 and blockNum == 3):
        sliceName = ['3','7','12','17','22','27','32','37']
    ADHip = '/AD_Hip' + str(blockNum)
    
    if blockNum == 1:
        ADHip = ''

    meanSliceNum = (np.sum(range(1,16)))/15.0 # based on z input found in align_v01_diffeo_by_katie_single.m
    numSlices = 15
    if (brainNum == 2 and blockNum == 3):
        meanSliceNum = (np.sum(range(1,8)))/7.0
        numSlices = 7
    if (brainNum == 3):
        meanSliceNum = (np.sum(range(1,14)))/13.0
        numSlices = 13
        if (blockNum == 2):
            meanSliceNum = (np.sum(range(1,13)))/12.0
            numSlices = 12
        elif (blockNum == 4):
            meanSliceNum = (np.sum(range(1,12)))/11.0
            numSlices = 11
    if (brainNum == 5):
        meanSliceNum = (np.sum(range(1,14)))/13.0
        numSlices = 13
        if (blockNum == 2):
            meanSliceNum = (np.sum(range(1,15)))/14.0
            numSlices = 14
        elif (blockNum == 3):
            meanSliceNum = (np.sum(range(1,9)))/8.0
            numSlices = 8

    mriParams = '/cis/home/kstouff4/Documents/HistoMRIPipeline/subject_Brain' + str(brainNum) + '_block_' + str(blockNum) + '_alignment_diffeo_test9_by_katie_single_round2/500_iters_' + dateNew + '/relTrans2.mat'
    mriTrans = '/cis/home/kstouff4/Documents/HistoMRIPipeline/subject_Brain'+ str(brainNum) +'_block_' + str(blockNum) + '_alignment_diffeo_test9_by_katie_single_round2/500_iters_' + dateNew + '/relTrans2.mat'
    mriCoords = '/cis/home/kstouff4/Documents/HistoMRIPipeline/subject_Brain'+ str(brainNum) +'_block_' + str(blockNum) + '_alignment_diffeo_test9_by_katie_single/coords.mat'
    
    if ('F' in dateNew):
        mriParams = '/cis/home/kstouff4/Documents/HistoMRIPipeline/subject_Brain' + str(brainNum) + '_block_' + str(blockNum) + '_alignment_diffeo_test9_by_katie_single/500_iters/relTrans.mat'
        mriTrans = '/cis/home/kstouff4/Documents/HistoMRIPipeline/subject_Brain'+ str(brainNum) +'_block_' + str(blockNum) + '_alignment_diffeo_test9_by_katie_single/500_iters/relTrans.mat'
        mriCoords = '/cis/home/kstouff4/Documents/HistoMRIPipeline/subject_Brain'+ str(brainNum) +'_block_' + str(blockNum) + '_alignment_diffeo_test9_by_katie_single/coords.mat'

    # Default to Tau
    if (brainNum == 2 and blockNum < 3):
        Ihisto = ['/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/histology/down_000/AD_Hip' + str(blockNum) + ADHip + '/Tau/BRC2614_AD_Block' + str(blockNum) + '_PHF-1_Location_5_corrected.tif',
                  '/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/histology/down_000/AD_Hip' + str(blockNum) + ADHip + '/Tau/BRC2614_AD_Block' + str(blockNum) + '_PHF-1_Location_6_corrected.tif',
                  '/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/histology/down_000/AD_Hip' + str(blockNum) + ADHip + '/Tau/BRC2614_AD_Block' + str(blockNum) + '_PHF-1_Location_7_corrected.tif',
                  '/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/histology/down_000/AD_Hip' + str(blockNum) + ADHip + '/Tau/BRC2614_AD_Block' + str(blockNum) + '_PHF-1_Location_8_corrected.tif',
                  '/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/histology/down_000/AD_Hip' + str(blockNum) + ADHip + '/Tau/BRC2614_AD_Block' + str(blockNum) + '_PHF-1_Location_9_corrected.tif',
                  '/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/histology/down_000/AD_Hip' + str(blockNum) + ADHip + '/Tau/BRC2614_AD_Block' + str(blockNum) + '_PHF-1_Location_10_corrected.tif',
                  '/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/histology/down_000/AD_Hip' + str(blockNum) + ADHip + '/Tau/BRC2614_AD_Block' + str(blockNum) + '_PHF-1_Location_11_corrected.tif',
                  '/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/histology/down_000/AD_Hip' + str(blockNum) + ADHip + '/Tau/BRC2614_AD_Block' + str(blockNum) + '_PHF-1_Location_12_corrected.tif',
                  '/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/histology/down_000/AD_Hip' + str(blockNum) + ADHip + '/Tau/BRC2614_AD_Block' + str(blockNum) + '_PHF-1_Location_13_corrected.tif',
                  '/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/histology/down_000/AD_Hip' + str(blockNum) + ADHip + '/Tau/BRC2614_AD_Block' + str(blockNum) + '_PHF-1_Location_14_corrected.tif',
                  '/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/histology/down_000/AD_Hip' + str(blockNum) + ADHip + '/Tau/BRC2614_AD_Block' + str(blockNum) + '_PHF-1_Location_15_corrected.tif',
                  '/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/histology/down_000/AD_Hip' + str(blockNum) + ADHip + '/Tau/BRC2614_AD_Block' + str(blockNum) + '_PHF-1_Location_16_corrected.tif',
                  '/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/histology/down_000/AD_Hip' + str(blockNum) + ADHip + '/Tau/BRC2614_AD_Block' + str(blockNum) + '_PHF-1_Location_17_corrected.tif',
                  '/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/histology/down_000/AD_Hip' + str(blockNum) + ADHip + '/Tau/BRC2614_AD_Block' + str(blockNum) + '_PHF-1_Location_18_corrected.tif',
                  '/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/histology/down_000/AD_Hip' + str(blockNum) + ADHip + '/Tau/BRC2614_AD_Block' + str(blockNum) + '_PHF-1_Location_19_corrected.tif'
        ]

        if (stainNum == 0):
            Ihisto = ['/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/histology/down_000/AD_Hip' + str(blockNum) + ADHip + '/Amyloid/BRC2614 AD Block' + str(blockNum) + ' 6E10 Location 5_corrected.tif',
                  '/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/histology/down_000/AD_Hip' + str(blockNum) + ADHip + '/Amyloid/BRC2614 AD Block' + str(blockNum) + ' 6E10 Location 6_corrected.tif',
                  '/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/histology/down_000/AD_Hip' + str(blockNum) + ADHip + '/Amyloid/BRC2614 AD Block' + str(blockNum) + ' 6E10 Location 7_corrected.tif',
                  '/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/histology/down_000/AD_Hip' + str(blockNum) + ADHip + '/Amyloid/BRC2614 AD Block' + str(blockNum) + ' 6E10 Location 8_corrected.tif',
                  '/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/histology/down_000/AD_Hip' + str(blockNum) + ADHip + '/Amyloid/BRC2614 AD Block' + str(blockNum) + ' 6E10 Location 9_corrected.tif',
                  '/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/histology/down_000/AD_Hip' + str(blockNum) + ADHip + '/Amyloid/BRC2614 AD Block' + str(blockNum) + ' 6E10 Location 10_corrected.tif',
                  '/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/histology/down_000/AD_Hip' + str(blockNum) + ADHip + '/Amyloid/BRC2614 AD Block' + str(blockNum) + ' 6E10 Location 11_corrected.tif',
                  '/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/histology/down_000/AD_Hip' + str(blockNum) + ADHip + '/Amyloid/BRC2614 AD Block' + str(blockNum) + ' 6E10 Location 12_corrected.tif',
                  '/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/histology/down_000/AD_Hip' + str(blockNum) + ADHip + '/Amyloid/BRC2614 AD Block' + str(blockNum) + ' 6E10 Location 13_corrected.tif',
                  '/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/histology/down_000/AD_Hip' + str(blockNum) + ADHip + '/Amyloid/BRC2614 AD Block' + str(blockNum) + ' 6E10 Location 14_corrected.tif',
                  '/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/histology/down_000/AD_Hip' + str(blockNum) + ADHip + '/Amyloid/BRC2614 AD Block' + str(blockNum) + ' 6E10 Location 15_corrected.tif',
                  '/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/histology/down_000/AD_Hip' + str(blockNum) + ADHip + '/Amyloid/BRC2614 AD Block' + str(blockNum) + ' 6E10 Location 16_corrected.tif',
                  '/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/histology/down_000/AD_Hip' + str(blockNum) + ADHip + '/Amyloid/BRC2614 AD Block' + str(blockNum) + ' 6E10 Location 17_corrected.tif',
                  '/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/histology/down_000/AD_Hip' + str(blockNum) + ADHip + '/Amyloid/BRC2614 AD Block' + str(blockNum) + ' 6E10 Location 18_corrected.tif',
                  '/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/histology/down_000/AD_Hip' + str(blockNum) + ADHip + '/Amyloid/BRC2614 AD Block' + str(blockNum) + ' 6E10 Location 19_corrected.tif'
            ]

    elif (brainNum == 2 and blockNum == 3):
        Ihisto = ['/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/histology/down_000/AD_Hip' + str(blockNum) + ADHip + '/Tau/BRC2614_AD_Block' + str(blockNum) + '_PHF-1_Location_5_corrected.tif',
                  '/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/histology/down_000/AD_Hip' + str(blockNum) + ADHip + '/Tau/BRC2614_AD_Block' + str(blockNum) + '_PHF-1_Location_6_corrected.tif',
                  '/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/histology/down_000/AD_Hip' + str(blockNum) + ADHip + '/Tau/BRC2614_AD_Block' + str(blockNum) + '_PHF-1_Location_7_corrected.tif',
                  '/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/histology/down_000/AD_Hip' + str(blockNum) + ADHip + '/Tau/BRC2614_AD_Block' + str(blockNum) + '_PHF-1_Location_8_corrected.tif',
                  '/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/histology/down_000/AD_Hip' + str(blockNum) + ADHip + '/Tau/BRC2614_AD_Block' + str(blockNum) + '_PHF-1_Location_9_corrected.tif',
                  '/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/histology/down_000/AD_Hip' + str(blockNum) + ADHip + '/Tau/BRC2614_AD_Block' + str(blockNum) + '_PHF-1_Location_10_corrected.tif',
                  '/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/histology/down_000/AD_Hip' + str(blockNum) + ADHip + '/Tau/BRC2614_AD_Block' + str(blockNum) + '_PHF-1_Location_11_corrected.tif'
        ]


        if (stainNum == 0):
            Ihisto = ['/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/histology/down_000/AD_Hip' + str(blockNum) + ADHip + '/Amyloid/BRC2614 AD Block' + str(blockNum) + ' 6E10 Location 5_corrected.tif',
                  '/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/histology/down_000/AD_Hip' + str(blockNum) + ADHip + '/Amyloid/BRC2614 AD Block' + str(blockNum) + ' 6E10 Location 6_corrected.tif',
                  '/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/histology/down_000/AD_Hip' + str(blockNum) + ADHip + '/Amyloid/BRC2614 AD Block' + str(blockNum) + ' 6E10 Location 7_corrected.tif',
                  '/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/histology/down_000/AD_Hip' + str(blockNum) + ADHip + '/Amyloid/BRC2614 AD Block' + str(blockNum) + ' 6E10 Location 8_corrected.tif',
                  '/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/histology/down_000/AD_Hip' + str(blockNum) + ADHip + '/Amyloid/BRC2614 AD Block' + str(blockNum) + ' 6E10 Location 9_corrected.tif',
                  '/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/histology/down_000/AD_Hip' + str(blockNum) + ADHip + '/Amyloid/BRC2614 AD Block' + str(blockNum) + ' 6E10 Location 10_corrected.tif',
                  '/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/histology/down_000/AD_Hip' + str(blockNum) + ADHip + '/Amyloid/BRC2614 AD Block' + str(blockNum) + ' 6E10 Location 11_corrected.tif'
            ]
    if (brainNum == 2):  
        histTrans = []
        for i in range(sliceNums[int(brainNum-1),int(blockNum-1),int(stainNum)]):
            if (stainNum == 0):
                histTrans.append('/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/histology/down_000/AD_Hip' + str(blockNum) + '/outputFigs/' + dateNew + '/Amyloid/Amyloid-BRC2614 AD Block' + str(blockNum) + ' 6E10 Location ' + str(i+5) + '_corrected.tif-to-L' + str(i+5) + mriSuffix + '.tif-vars.npz')
            elif (stainNum == 2):
                histTrans.append('/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/histology/down_000/AD_Hip' + str(blockNum) + '/outputFigs/' + dateNew + '/Tau/Tau-BRC2614_AD_Block' + str(blockNum) + '_PHF-1_Location_' + str(i+5) + '_corrected.tif-to-L' + str(i+5) + mriSuffix + '.tif-vars.npz')

    else:
        Ihisto = []
        histTrans = []
        for i in range(sliceNums[int(brainNum-1),int(blockNum-1),int(stainNum)]):
            Ihisto.append('/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain' + str(brainNum) + '/histology/down_000/AD_Hip' + str(blockNum) + '/' + stainString + '/Brain ' + str(brainNum) + '-Block ' + str(blockNum) + ' L' + sliceName[i] + ' ' + stain[stainNum] + '_crop_filterNS.tif')
            histTrans.append('/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain' + str(brainNum) + '/histology/down_000/AD_Hip' + str(blockNum) + '/outputFigs/' + dateNew + '/' + stainString + '/' + stainString + '-Brain ' + str(brainNum) + '-Block ' + str(blockNum) + ' L' + sliceName[i] + ' ' + stain[stainNum] + '_crop_filterNS.tif-to-L' + sliceName[i] + mriSuffix + '.tif-vars.npz')

    #x0m,x1m,mriI = load_image_only('/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain5/histology/alignment/block2/itk/L10.tif',0.06656,0.06656)        
    #X0mP,X1mP = np.meshgrid(x0m,x1m,indexing='ij')
    X0mP = None
    X1mP = None

    runPipelineForBlock(brainNum, blockNum-1, stainNum, meanSliceNum, Ihisto, mriParams, histTrans, mriTrans, mriCoords,X0mP,X1mP,saveName='_'+mriSuffix+'_' + dateNew,lowBound=-1,highBound=20,numSlices=numSlices)
    return

def getCompositionManual(brainNum,blockNum,sliceNum,stain,dateNew,mriSuffix,matFile):
    
    if ('XM' in matFile):
        XMstring = 'XM'
    else:
        XMstring = ''
    print('in composition manual')
    stainNum = 2 # 0 indexed: Amyloid, KFB, Tau
    stainString = "Tau"
    if (stain == 'Amyloid'):
        stainString = "Amyloid"
        stainNum = 0

    stain = ['6E10','LFB H&E','PHF-1']
    sliceName = ['1','5','10','15','20','25','30','35','40','45','50','55','60','65','70','75']
    
    if (brainNum == 5 and blockNum == 3):
        sliceName = ['3','7','12','17','22','27','32','37']
    Ihisto = []
    histTrans = []
        
    if (brainNum == 2):
        ADHip = '/AD_Hip' + str(blockNum)
        if (blockNum == 1):
            ADHip = ''
            end = 15
        if (blockNum == 2):
            end = 15
        if (blockNum == 3):
            end = 7
        for sl in range(end):
            if (stainNum == 0):
                Ihisto.append('/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/histology/down_000/AD_Hip' + str(blockNum) + ADHip + '/Amyloid/BRC2614 AD Block' + str(blockNum) + ' 6E10 Location ' + str(sl + 5) + '_corrected.tif')
                histTrans.append('/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/histology/down_000/AD_Hip' + str(blockNum) + '/outputFigs/' + dateNew + '/Amyloid/Amyloid-BRC2614 AD Block' + str(blockNum) + ' 6E10 Location ' + str(sl+5) + '_corrected.tif-to-L' + str(sl+5) + mriSuffix + '.tif-vars_manual' + XMstring + '.npz')
            elif (stainNum == 2):
                Ihisto.append('/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/histology/down_000/AD_Hip' + str(blockNum) + ADHip + '/Tau/BRC2614_AD_Block' + str(blockNum) + '_PHF-1_Location_' + str(sl+5) + '_corrected.tif')
                histTrans.append('/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/histology/down_000/AD_Hip' + str(blockNum) + '/outputFigs/' + dateNew + '/Tau/Tau-BRC2614_AD_Block' + str(blockNum) + '_PHF-1_Location_' + str(sl+5) + '_corrected.tif-to-L' + str(sl+5) + mriSuffix + '.tif-vars_manual' + XMstring + '.npz')
         
        
    else:
        if (blockNum == 1):
            end = 13
        if (blockNum == 2):
            end = 14
        if (blockNum == 3):
            end = 8
        for sl in range(end):
            Ihisto.append('/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain' + str(brainNum) + '/histology/down_000/AD_Hip' + str(blockNum) + '/' + stainString + '/Brain ' + str(brainNum) + '-Block ' + str(blockNum) + ' L' + sliceName[sl] + ' ' + stain[stainNum] + '_crop_filterNS.tif')
            histTrans.append('/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain' + str(brainNum) + '/histology/down_000/AD_Hip' + str(blockNum) + '/outputFigs/' + dateNew + '/' + stainString + '/' + stainString + '-Brain ' + str(brainNum) + '-Block ' + str(blockNum) + ' L' + sliceName[sl] + ' ' + stain[stainNum] + '_crop_filterNS.tif-to-L' + sliceName[sl] + mriSuffix + '.tif-vars_manual' + XMstring + '.npz')
    
    meanSliceNum = (np.sum(range(1,16)))/15.0 # based on z input found in align_v01_diffeo_by_katie_single.m
    numSlices = 15
    if (brainNum == 2 and blockNum == 3):
        meanSliceNum = (np.sum(range(1,8)))/7.0
        numSlices = 7
    if (brainNum == 3):
        meanSliceNum = (np.sum(range(1,14)))/13.0
        numSlices = 13
        if (blockNum == 2):
            meanSliceNum = (np.sum(range(1,13)))/12.0
            numSlices = 12
        elif (blockNum == 4):
            meanSliceNum = (np.sum(range(1,12)))/11.0
            numSlices = 11
    if (brainNum == 5):
        meanSliceNum = (np.sum(range(1,14)))/13.0
        numSlices = 13
        if (blockNum == 2):
            meanSliceNum = (np.sum(range(1,15)))/14.0
            numSlices = 14
        elif (blockNum == 3):
            meanSliceNum = (np.sum(range(1,9)))/8.0
            numSlices = 8

    mriParams = '/cis/home/kstouff4/Documents/HistoMRIPipeline/subject_Brain' + str(brainNum) + '_block_' + str(blockNum) + '_alignment_diffeo_test9_by_katie_single_round2/500_iters_' + dateNew + '/relTrans2.mat'
    mriTrans = '/cis/home/kstouff4/Documents/HistoMRIPipeline/subject_Brain'+ str(brainNum) +'_block_' + str(blockNum) + '_alignment_diffeo_test9_by_katie_single_round2/500_iters_' + dateNew + '/relTrans2.mat'
    mriCoords = '/cis/home/kstouff4/Documents/HistoMRIPipeline/subject_Brain'+ str(brainNum) +'_block_' + str(blockNum) + '_alignment_diffeo_test9_by_katie_single/coords.mat'
    
    if ('F' in dateNew):
        mriParams = '/cis/home/kstouff4/Documents/HistoMRIPipeline/subject_Brain' + str(brainNum) + '_block_' + str(blockNum) + '_alignment_diffeo_test9_by_katie_single/500_iters/relTrans.mat'
        mriTrans = '/cis/home/kstouff4/Documents/HistoMRIPipeline/subject_Brain'+ str(brainNum) +'_block_' + str(blockNum) + '_alignment_diffeo_test9_by_katie_single/500_iters/relTrans.mat'
        mriCoords = '/cis/home/kstouff4/Documents/HistoMRIPipeline/subject_Brain'+ str(brainNum) +'_block_' + str(blockNum) + '_alignment_diffeo_test9_by_katie_single/coords.mat'

    print(matFile)
    if (path.exists(matFile)):
        if ('XM' in matFile):
            print('XM in matfile')
            if ('npz' in matFile):
                params = np.load(matFile)
                X0mP = params['X0m']
                X1mP = params['X1m']
                phi = np.stack((X0mP,X1mP),axis=-1)
                X0I = X0mP
                X1I = X0mP
                A = np.eye(3)
                np.savez(histTrans[sliceNum],A=A,phi=phi,X0I=X0I,X1I=X1I,X0m=X0mP,X1m=X1mP)
            else:
                params=sp.io.loadmat(matFile,appendmat=False,struct_as_record=False)
                X0mP = params['X0M'].astype('float32')
                X1mP = params['X1M'].astype('float32')
        elif ('npz' in matFile):
            params = np.load(matFile)
            phi = params['phi']
            A = params['A']
            X0I = params['X0I']
            X1I = params['X1I']
            np.savez(histTrans[sliceNum],A=A,phi=phi,X0I=X0I,X1I=X1I)
            X0mP = None
            X1mP = None
        else:
            params = sp.io.loadmat(matFile,appendmat=False, struct_as_record=False)
            A = params['B'].astype('float32')
            Ac = np.copy(A)
            Ac[0,1] = A[1,0]
            Ac[1,0] = A[0,1]
            Ac[0,2] = A[1,2]
            Ac[1,2] = A[0,2]
            #A = np.linalg.inv(Ac)
            A = Ac
            phi = params['phi'].astype('float32')
            print('phi shape should be R x C by 2')
            print(phi.shape)
            phic = np.copy(phi)
            phic[...,0] = phi[...,1]
            phic[...,1] = phi[...,0]
            phi = phic
            X0I = params['YJ'].astype('float32')
            X1I = params['XJ'].astype('float32')
            print('min of coordinates')
            print(np.min(X0I))
            print(np.min(X1I))
            np.savez(histTrans[sliceNum],A=A,phi=phi,X0I=X0I,X1I=X1I)
            #X0mP = phi[...,0]
            #X1mP = phi[...,1]
            X0mP = None
            X1mP = None
    else:
        X0mP = None
        X1mP = None

    runPipelineForBlock(brainNum, blockNum-1, stainNum, meanSliceNum, Ihisto, mriParams, histTrans, mriTrans, mriCoords,X0mP,X1mP,saveName='_'+mriSuffix+'_' + dateNew,lowBound=sliceNum-1,highBound=sliceNum+1,numSlices=numSlices)
    return