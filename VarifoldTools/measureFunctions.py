################################################################
# Author: Katie Stouffer, 2/6/22

# Tools for constructing, manipulating, sampling, etc. varifolds. These tools come mostly from the General_Varifold and Varifold framework python notebooks.

################################################################

import sys
sys.path.append('../')
sys.path.append('../SurfaceTools/')

import deformSegmentations3D as ds
import vtkFunctions as vt
import varifoldFunctions as vf

# Libraries to Import
import numpy as np
import math
import sys
import matplotlib.pyplot as plt
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

import seaborn as sns
import pandas as pd
import random

###############################################################
# Functions for Generating Measure Representation

def crossProduct(a,b):
    p0 = a[...,1]*b[...,2] - a[...,2]*b[...,1]
    p1 = a[...,2]*b[...,0] - a[...,0]*b[...,2]
    p2 = a[...,0]*b[...,1] - a[...,1]*b[...,0]
    
    return p0,p1,p2

def dotProduct(a,b):
    p0 = a[...,0]*b[...,0] + a[...,1]*b[...,1] + a[...,2]*b[...,2]
    return p0
    
def getAreaWeightsExvivo(brainNum,blocks,stainNum,mriSuffix,dateNew):
    '''
    Calculates weights of particles (pixels) in 3D space as avg of areas of 4 polygons that a given point is vertex to
    in theory, \sum w_j \delta(xj) \cross \delta(f_xj) --> \sum w_j |d\phi(xj)| \delta(phi(xj)) \cross \delta(f_\phi(xj))
    here, \phi is NOT a diffeomorphism strictly speaking because it is not invertible, however it represents how much a 
    unit cube changes under the transformation phi. 
    
    Here, w_j is defined as dxdy = (1/4)(dxdy + dxdy + dxdy + dxdy) as we assume uniformity in pixel size initially
    Under the trasnformation, we assume the pixels are close enough so that they can still be considered a 2D manifold (locally similar to R^2)
    
    Effectively, we consider these particles to still be 2D in that \int_x\int_y delta(X) dxdy = 1, or more simply \int_x \delta(X) dx = 1 (considering x multidimensional)
    , thus, under teh composition of transformations, we have \int_(phi(x)) delta(phi(X)) |dphi| dx and tehrefore, 
    we can use teh same formula to think of w_j as |dphi|dx = "area" of particle
    '''
    savename = '_'+mriSuffix+'_' + dateNew
    for blockNum in blocks:
        high = 15
        if (brainNum == 2 and blockNum == 2):
            high = 7
        for sl in range(0,high):
            locationNum = sl+5
            f = '/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/2DHisto_to_MRI_coords_block' + str(blockNum+1) + '_stain' + str(stainNum) + '_slice' + str(sl) + savename + '.npz'
            if (brainNum == 5):
                f = '/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/2DHisto_to_MRI_coords_block' + str(blockNum+1) + '_stain' + str(stainNum) + '_slice' + str(sl) + savename + '.npz'
            if (brainNum == 3):
                f = '/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/2DHisto_to_MRI_coords_block' + str(blockNum+1) + '_stain' + str(stainNum) + '_slice' + str(sl) + savename + '.npz'
            if (not path.exists(f)):
                continue
            histocoords_3D_npz = np.load(f) #_newmai.npz for brain2
            histocoords_3D = histocoords_3D_npz['coordsT_mai']
            histocoords_3D = np.squeeze(histocoords_3D).astype('float32')
            histocoords_3D_npz.close()
        
            W = np.zeros((histocoords_3D.shape[0],histocoords_3D.shape[1]))
            h0 = np.zeros((histocoords_3D.shape[0]+1,histocoords_3D.shape[1]+1,histocoords_3D.shape[2]))
            h1 = np.zeros_like(h0)
            h2 = np.zeros_like(h0)
            h3 = np.zeros_like(h0)
            
            h0[1:,1:,:] = histocoords_3D # shifted down 1 and to right 1
            h1[1:,0:-1,:] = histocoords_3D # shifted left 1
            h2[0:-1,0:-1,:]=histocoords_3D
            h3[0:-1,1:,:] = histocoords_3D # shifted up 1
            '''
            # old
            a0,a1,a2 = crossProduct(h1-h0,h2-h0)
            b0,b1,b2 = crossProduct(h2-h0,h3-h0)
            c0,c1,c2 = crossProduct(h3-h0,h1-h0)
            
            # get areas 
            s0 = a0+b0+c0
            s1 = a1+b1+c1
            s2 = a2+b2+c2
            
            A = 0.5*np.sqrt(s0**2 + s1**2 + s2**2) # 1/2 magnitude of normal
            '''
            
            a0,a1,a2 = crossProduct(h1-h0,h2-h0)
            b0,b1,b2 = crossProduct(h2-h0,h3-h0)
            tri1 = a0**2 + a1**2 + a2**2
            tri1 = np.sqrt(tri1)/2.0
            tri2 = b0**2 + b1**2 + b2**2
            tri2 = np.sqrt(tri2)/2.0
            A = tri1+tri2
            W = A[1:,1:] # should be area except for rightmost and bottommost pixels, but should be area
            np.savez('/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/2DHisto_to_MRI_coords_block' + str(blockNum+1) + '_stain' + str(stainNum) + '_slice' + str(sl) + savename + '_weights2.npz',W=W)
            
    return

def downsampleExVivoAndConvert(brainNum,blocks,numLab,stainNum,exp,mriSuffix,dateNew,redone,byLabel=True,histoLab=False,weightBy=2,start=0,end=15):
    '''
    Downsamples fine particle locations to coarse particle locations using intrinsic "Haar" kernel where each 
    fine particle is weighted with 1 or 0 per coarse particle.
    
    Downsample by 2^exp
    
    Saves format for both varifold and general varifold 
    
    if brainNum != 2, use background and crop mask to ignore regions where tissue was damaged 
    
    mriSuffix and dateNew both strings (e.g. '_pil' and '020222'
    redone = redone_11222
    '''

    savename = '_'+mriSuffix+'_' + dateNew
    if (stainNum == 0):
        stainString = "Amyloid"
        stainRef = "6E10"
    elif (stainNum == 2):
        stainString = "Tau"
        stainRef = "PHF-1"
    
    if (brainNum == 2):
        sliceNames = [5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]
    elif (brainNum == 5):
        sliceNames = [1,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75]
        
    # For checking conservation of mass
    weightTotal = 0 # mm^2
    tauTotal = 0 # numParticles
    labTotal = np.zeros((numLab,1)) # mm^2 per region
    labTotalTau = np.zeros((numLab,1))
    weightTotalA = 0
    tauTotalA = 0
    labTotalA = np.zeros((numLab,1))
    labTotalTauA = np.zeros((numLab,1))
    finePartCoords = [] # x,y,z coords
    finePartTau = [] # tau particles in tissue in fine particle (0 or 1)
    finePartLab = [] # vector valued fraction of area in region (0 or 1)
    finePartTauLab = [] # vector valued tau in region (0 or 1)
    finePartWeights = [] # 0 or 0.00208 mm^2 (area of pixel)
    finePartTauArea = []
    finePartTauOrient = []
    finePartTauRound = []
    
    for blockNum in blocks:
        if (brainNum == 2 and blockNum == 1):
            start = 1
            start = 0 # use all of the slices (03/04/22)
        if (brainNum == 2 and blockNum == 2):
            end = 7
            start = 0
        if (brainNum == 5 and blockNum == 2):
            end = 8
            sliceNames = [3,7,12,17,22,27,32,37]
        if (brainNum == 3):
            end = 12
            if (blockNum == 1):
                end = 13
            elif (blockNum == 3):
                end = 11
        for sl in range(start,end):
            locationNum = sliceNames[sl]
            fileHi = '/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/2DHisto_to_MRI_coords_block' + str(blockNum+1) + '_stain' + str(stainNum) + '_slice' + str(sl) + savename + '.npz'
            if (brainNum == 5):
                fileHi = '/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/2DHisto_to_MRI_coords_block' + str(blockNum+1) + '_stain' + str(stainNum) + '_slice' + str(sl) + savename + '.npz'
            if (brainNum == 3):
                fileHi = '/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/2DHisto_to_MRI_coords_block' + str(blockNum+1) + '_stain' + str(stainNum) + '_slice' + str(sl) + savename + '.npz'
            if (not path.exists(fileHi)):
                continue;
            histocoords_3D_npz = np.load(fileHi)
            if (brainNum == 2):
                locationNum = sl+5
            histocoords_3D = histocoords_3D_npz['coordsT_mai']
            if (weightBy < 2):
                histocoords_3D = histocoords_3D_npz['coordsH'] # assume interested in histology space
            histocoords_3D = histocoords_3D.astype('float32')
            histocoords_3D_npz.close()
            if (weightBy == 2):
                wF = np.load('/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/2DHisto_to_MRI_coords_block' + str(blockNum+1) + '_stain' + str(stainNum) + '_slice' + str(sl) + savename + '_weights2.npz')
                weightsFine = wF['W'].astype('float32')
            elif (weightBy == 3):
                wF = np.load('/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/2DHisto_to_MRI_coords_block' + str(blockNum+1) + '_stain' + str(stainNum) + '_slice' + str(sl) + savename + '_weights3D.npz')
                weightsFine = wF['W'].astype('float32')
            else:
                weightsFine = np.ones((histocoords_3D.shape[0],histocoords_3D.shape[1]))*0.00208*0.00208 # default to area only
            #weightsFine = np.ones((histocoords_3D.shape[0],histocoords_3D.shape[1]))
            
            print("sum of weights fine before " + str(np.sum(weightsFine)))
            backGroundWeights = np.load('/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/Block' + str(blockNum + 1) + '/'+stainString + '/BackgroundMasks/down_000/Location_' + str(locationNum) + '.npz')
            backGroundWeights = backGroundWeights['openmaskup']
            backGroundWeights[0:50,...] = 1.0
            backGroundWeights[:,0:50] = 1.0
            backGroundWeights[-50:,:] = 1.0
            backGroundWeights[:,-50:] = 1.0
            backGroundWeights = (backGroundWeights == 0).astype('float32')
            if (brainNum != 2):
                print("reading in mask")
                maskWeights = np.load('/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain' + str(brainNum) + '/histology/down_000/AD_Hip' + str(blockNum + 1) + '/'+stainString + '/Brain '+ str(brainNum) + '-Block ' + str(blockNum+1) + ' L' + str(locationNum) + ' ' + stainRef + '_cropMask.npz')
                maskWeights = maskWeights['mask'].astype('float32') # 1 for black pixels
                backGroundWeights = backGroundWeights * (1.0 - maskWeights) # zero out black pixels
                print(np.unique(maskWeights))
            weightsFine = weightsFine * backGroundWeights # zero out background (+/- black pixels)
            f,ax = plt.subplots()
            im = ax.imshow(weightsFine>0,origin = 'lower', cmap='viridis')
            f.colorbar(im,ax=ax)
            ax.set_title("Before downsampling")
            f.canvas.draw()
            if (brainNum == 2):
                if (stainNum == 0):
                        descriptionFile = np.genfromtxt('/cis/home/dtward/Public/For_Katie/From_Paige/brain2_block' + str(blockNum+1) + '_location' + str(sl+5) + '/description.csv',delimiter=',')
                elif (stainNum == 2):
                        #descriptionFile = np.genfromtxt('/cis/home/slee508/my_documents/Brain' + str(brainNum) + '/asset/block' + str(blockNum+1) + '/description_' + str(sl+5) + '.csv',delimiter=',')
                        descriptionFile = np.genfromtxt('/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/ModelOutput/Brain' + str(brainNum) + '/AD_Hip' + str(blockNum+1) + '/Tau/Location_' + str(sl+5) + '_description.csv',delimiter=',')
            elif (brainNum == 5):
                if (stainNum == 0):
                        descriptionFile = np.genfromtxt('/cis/home/dtward/Public/For_Katie/From_Paige/brain2_block' + str(blockNum+1) + '_location' + str(sl+5) + '/description.csv',delimiter=',')
                elif (stainNum == 2):
                        #descriptionFile = np.genfromtxt('/cis/home/slee508/my_documents/Brain' + str(brainNum) + '/asset/block' + str(blockNum+1) + '/description_' + str(sl+5) + '.csv',delimiter=',')
                        descriptionFile = np.genfromtxt('/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/ModelOutput/Brain' + str(brainNum) + '/AD_Hip' + str(blockNum+1) + '/Tau/L' + str(locationNum) + '_description.csv',delimiter=',')
            elif (brainNum == 3):
                if (stainNum == 2):
                    descriptionFile = np.genfromtxt('/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/ModelOutput/Brain' + str(brainNum) + '/AD_Hip' + str(blockNum+1) + '/Tau/L' + str(locationNum) + '_description.csv',delimiter=',')
            else:
                if (stainNum == 2):
                    if (brainNum == 5):
                        if ((blockNum == 0 and locationNum == 15) or (blockNum == 1 and (locationNum == 10 or locationNum == 15))):
                            descriptionFile = np.genfromtxt('/cis/home/slee508/my_documents/Brain' + str(brainNum) + '/asset/block' + str(blockNum+1) + '/pmap_Brain ' + str(brainNum) + '-Block ' + str(blockNum+1) + ' L' + str(locationNum) + ' PHF-1_crop.csv',delimiter=',')
                        else:
                            descriptionFile = np.genfromtxt('/cis/home/slee508/my_documents/Brain' + str(brainNum) + '/asset/block' + str(blockNum+1) + '/description_pmap_Brain ' + str(brainNum) + '-Block ' + str(blockNum+1) + ' L' + str(locationNum) + ' PHF-1_crop.csv',delimiter=',')
                    else:
                        descriptionFile = np.genfromtxt('/cis/home/slee508/my_documents/Brain' + str(brainNum) + '/asset/block' + str(blockNum+1) + '/description_pmap_Brain ' + str(brainNum) + '-Block ' + str(blockNum+1) + ' L' + str(locationNum) + ' PHF-1_crop.csv',delimiter=',')

            histoTau = np.zeros((histocoords_3D.shape[0],histocoords_3D.shape[1]))
            histoTauArea = np.zeros((histocoords_3D.shape[0],histocoords_3D.shape[1]))
            histoTauOrient = np.zeros((histocoords_3D.shape[0],histocoords_3D.shape[1]))
            histoTauRound = np.zeros((histocoords_3D.shape[0],histocoords_3D.shape[1]))
            
            rows = descriptionFile.shape[0]
            # Regenerate graphs with labels 
            if (brainNum == 2):
                labelIm_vars = np.load('/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/2DHisto_to_MRI_coords_block' + str(blockNum+1) + '_stain' + str(stainNum) + '_loc' + str(locationNum) + '_Newlabels_redone_11222' + savename + '.npz')
            if (brainNum == 5):
                labelIm_vars = np.load('/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/2DHisto_to_MRI_coords_block' + str(blockNum+1) + '_stain' + str(stainNum) + '_loc' + str(locationNum) + '_Newlabels_' + redone + savename + '.npz')
            elif (brainNum == 3):
                labelIm_vars = np.load('/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/2DHisto_to_MRI_coords_block' + str(blockNum+1) + '_stain' + str(stainNum) + '_loc' + str(locationNum) + '_Newlabels_' + redone + savename + '.npz')
            else:
                labelIm_vars = np.load('/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/2DHisto_to_MRI_coords_block' + str(blockNum+1) + '_stain' + str(stainNum) + '_loc' + str(locationNum) + '_Newlabels_' + redone + savename + '.npz')
            labelIm = labelIm_vars['L_histo_Mode']
            labelImbin = multiToBin(labelIm,numLab) # translate into one channel per label (only 0s and 1s)
            labelIm_vars.close()
            labelImTau = np.zeros_like(labelImbin)
            
            print("sizes of binary label im and histotau")
            print(labelImbin.shape)
            print(labelImTau.shape)
            print(histoTau.shape)
            print(weightsFine.shape)
            print(histocoords_3D.shape)

            for r in range(1,rows):
                gCenterX = int(descriptionFile[r,7]) #(col coord??)
                gCenterY = int(descriptionFile[r,8])
                histoTau[gCenterY,gCenterX] = 1.0
                histoTauArea[gCenterY,gCenterX] = descriptionFile[r,9]
                histoTauOrient[gCenterY,gCenterX] = descriptionFile[r,10]
                histoTauRound[gCenterY,gCenterX] = descriptionFile[r,11]
                #weightsFine = 0.002**2*weightsFine # non-zero weight for all sampled TISSUE region
                # downsample histoTau and histocoords_3D
                  # Transition between level 2 and level 1 (intermediate)
            for lab in range(labelImbin.shape[-1]):
                labelImTau[...,lab] = histoTau[...,None]*labelImbin[...,lab] # only count ones that are in label and tau
            print("after added histoTau, should be equal")
            print(np.sum(histoTau))
            print(np.sum(labelImTau))
            
            '''
            Coarse particle coordinates = average
            labelIm = sum only
            histoTau = sum weightsFine*histoTau (assuming weightsFine = 0 or 1)
            weightsFine = sum only (counts number of tissue particles, to be multiplied in the end by 0.00208^2)
            '''
            #histoTauWeight = histoTau*weightsFine # zero out background particles
            tauTotal += np.sum(histoTau*backGroundWeights) # just counting non background particles
            if (True):
                #weightsFine = weightsFine * 0.00208**2
                histoTauWeight = histoTau*backGroundWeights # zero out background particles, function value is now particles*area
                histoTauAreaWeight = histoTauArea*backGroundWeights
                histoTauOrientWeight = histoTauOrient*backGroundWeights
                histoTauRoundWeight = histoTauRound*backGroundWeights
                weightTotal += np.sum(weightsFine)
                for lab in range(labelImbin.shape[-1]):
                    print("lab covered " + str(lab))
                    labelImbin[...,lab] = labelImbin[...,lab]*weightsFine[...,None] # zero out labels not from tissue sampled, function value is area of label sampled
                    labelImTau[...,lab] = labelImTau[...,lab]*backGroundWeights[...,None] # zero out labels not from tissue sampled, assume function value # tau particles
                    labTotal[lab] += np.sum(labelImbin[...,lab]) # number of pixels*area per pixel in this label
                    labTotalTau[lab] += np.sum(labelImTau[...,lab])
                    '''
                    wFd = 1./weightsFine
                    wFd[wFd == np.inf] = 0
                    labTotalTau[lab] += np.sum(labelImTau[...,lab]*wFd[...,None])
                    '''
                # This implementation seems to doublecount things, instead, try to just take all at once
                print("Values before")
                print("Tau total and contained in label Im")
                print(np.sum(histoTauWeight))
                print(np.sum(labelImTau))
                print("Weights total and contained in labelIm")
                print(np.sum(weightsFine))
                print(np.sum(labelImbin))
                print(np.sum(labelImbin,axis=(0,1,2)))
                print(np.sum(labelImTau,axis=(0,1,2)))
                down = 2**exp
                nx = histocoords_3D.shape[0]
                #Katie: floor division (discards fractional part)
                nxd = nx//down
                coordTemp = np.zeros((nxd,histocoords_3D.shape[1],histocoords_3D.shape[2],histocoords_3D.shape[3])) # assume 3D
                tauTemp = np.zeros((nxd,histoTauWeight.shape[1])) # assume 3D
                tauAreaTemp = np.zeros((nxd,histoTauWeight.shape[1])) # assume 3D
                tauOrientTemp = np.zeros((nxd,histoTauWeight.shape[1])) # assume 3D
                tauRoundTemp = np.zeros((nxd,histoTauWeight.shape[1])) # assume 3D
               
                labelImbinTemp = np.zeros((nxd,labelImbin.shape[1],labelImbin.shape[2],labelImbin.shape[3]))
                labelImTauTemp = np.zeros((nxd,labelImTau.shape[1],labelImTau.shape[2],labelImTau.shape[3]))
                weightsTemp = np.zeros((nxd,weightsFine.shape[1]))
                for i in range(down):
                    coordTemp += histocoords_3D[i:nxd*down:down,...]/down
                    tauTemp += histoTauWeight[i:nxd*down:down,...] # just add, don't divide, total number of tau particles
                    labelImbinTemp += labelImbin[i:nxd*down:down,...] # add, don't divide, total area in label
                    labelImTauTemp += labelImTau[i:nxd*down:down,...] # total tau particles in label
                    weightsTemp += weightsFine[i:nxd*down:down,...] # total area 
                    tauAreaTemp += histoTauAreaWeight[i:nxd*down:down,...] # total area of all tau tangles
                    tauOrientTemp += histoTauOrientWeight[i:nxd*down:down,...] # total orientation of all tau tangles
                    tauRoundTemp += histoTauRoundWeight[i:nxd*down:down,...] # total roundness of all tau tangles
                    
                nx = histocoords_3D.shape[1]
                nxd = nx//down
                coordTemp2 = np.zeros((coordTemp.shape[0],nxd,coordTemp.shape[2],coordTemp.shape[3])) # assume 3D
                tauTemp2 = np.zeros((tauTemp.shape[0],nxd)) # assume 3D
                labelImbinTemp2 = np.zeros((labelImbinTemp.shape[0],nxd,labelImbinTemp.shape[2],labelImbinTemp.shape[3]))
                labelImTauTemp2 = np.zeros((labelImTauTemp.shape[0],nxd,labelImTauTemp.shape[2],labelImTauTemp.shape[3]))
                weightsTemp2 = np.zeros((weightsTemp.shape[0],nxd))
                tauAreaTemp2 = np.zeros((tauAreaTemp.shape[0],nxd))
                tauOrientTemp2 = np.zeros((tauOrientTemp.shape[0],nxd))
                tauRoundTemp2 = np.zeros((tauRoundTemp.shape[0],nxd))
                
                for i in range(down):
                    coordTemp2 += coordTemp[:,i:nxd*down:down,...]/down
                    tauTemp2 += tauTemp[:,i:nxd*down:down,...] # just add, don't divide
                    tauAreaTemp2 += tauAreaTemp[:,i:nxd*down:down,...]
                    tauOrientTemp2 += tauOrientTemp[:,i:nxd*down:down,...]
                    tauRoundTemp2 += tauRoundTemp[:,i:nxd*down:down,...]
                    labelImbinTemp2 += labelImbinTemp[:,i:nxd*down:down,...] # add, don't divide
                    labelImTauTemp2 += labelImTauTemp[:,i:nxd*down:down,...]
                    weightsTemp2 += weightsTemp[:,i:nxd*down:down,...]        
  
                histocoords_3D = coordTemp2
                labelImbin = labelImbinTemp2
                labelImTau = labelImTauTemp2
                histoTauWeight = tauTemp2
                weightsFine = weightsTemp2
                
                histoTauAreaWeight = tauAreaTemp2 #/histoTauWeight # avg area (save total area for now)
                histoTauOrientWeight = tauOrientTemp2 #/histoTauWeight # avg orientation (save total )
                histoTauRoundWeight = tauRoundTemp2 #/histoTauWeight # avg roundness (save total)
                
                # do the area and orientation and roundness per tau tangles (so per weight and then per tau tangles)
                print("Values After")
                print(np.sum(histoTauWeight))
                print(np.sum(labelImTau))
                print("Weights total and contained in labelIm")
                print(np.sum(weightsFine))
                print(np.sum(labelImbin))
                print(np.sum(labelImbin,axis=(0,1,2)))
                print(np.sum(labelImTau,axis=(0,1,2)))

                # Normalize now because done with this coarse particle grid 
                # KATIE LEFT OFF HERE
                weightTotalA += np.sum(weightsFine)
                normWeights = 1./weightsFine # weightsFine now is .00208^2 * numpixels area
                normWeights[normWeights==np.inf] = 0.0
                histoTau = normWeights*histoTauWeight # htWght = numpixels=particles / * numpixels area particles per area sampled 
                histoTauArea = normWeights*histoTauAreaWeight
                histoTauOrient = normWeights*histoTauOrientWeight
                histoTauRound = normWeights*histoTauRoundWeight 
                tauTotalA += np.sum(histoTau*weightsFine)
                labelImTauSum = np.sum(labelImTau,axis=-1) # total tau in all labels
                print("should be equal (tau total) ")
                print(np.sum(labelImTauSum - histoTau[...,None]*weightsFine[...,None]))
                f,ax = plt.subplots()
                im = ax.imshow(np.squeeze(np.sum(labelImbin,axis=-1))-weightsFine)
                f.colorbar(im,ax=ax)
                f.canvas.draw()
                if (histoLab):
                    labelImTauSum = np.sum(labelImTau,axis=-1) # total tau in all labels
                    print("should be equal (tau total) ")
                    print(np.sum(labelImTauSum - histoTau*weightsFine))
                    labelIm_vars = np.load('/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/Block' + str(blockNum+1) + '/' + stainString + '/eXuSegmentations/Compiled/Location_' + str(locationNum) + '.npz')
                    labelImExu = labelIm_vars['labelImEx']
                    # Make label mask and tau and coordinates the same size (default to downsample size)
                    labelImExuResize = np.zeros_like(histoTau)
                    off1 = histocoords_3D.shape[0]-labelImExu.shape[0]
                    off2 = histocoords_3D.shape[1]-labelImExu.shape[1]
                    print("offsets are " + str(off1) + ", " + str(off2))
                    if (off1 <= 0):
                        if (off2 < 0):
                            labelImExuResize = labelImExu[0:histocoords_3D.shape[0],-off2:]
                        else:
                            labelImExuResize[:,off2:] = labelImExu[0:histocoords_3D.shape[0],:]
                    else:
                        if (off2 < 0):
                            labelImExuResize[0:-off1,:] = labelImExu[:,-off2:]
                        else:
                            labelImExuResize[0:-off1,off2:] = labelImExu
                    labelImExu = labelImExuResize[...,None]
                    labelImbin = multiToBin(labelImExu,numLab) # translate into one channel per label
                    labelIm_vars.close()
                    labelImTau = np.zeros_like(labelImbin)
                    for lab in range(numLab):
                        labelImTau[...,lab] = labelImbin[...,lab]*labelImTauSum*normWeights[...,None] # got rid of 1 / 0.00208^2 factor
                else:
                    for lab in range(numLab):
                        print("lab covered " + str(lab))
                        labTotalA[lab] += np.sum(labelImbin[...,lab])
                        labelImbin[...,lab] = labelImbin[...,lab]*normWeights[...,None] # now gives total amount of area in label per coarse particle / total area of coarse particle (to give fraction of area sampled in region)
                        labelImTau[...,lab] = labelImTau[...,lab]*normWeights[...,None] # total tau in label / area
                        normWeights2 = 1./labelImbin[...,lab]
                        normWeights2[normWeights2==np.inf] = 0.0
                        #labelImTau[...,lab] = labelImTau[...,lab]*(1/(0.00208**2)) # TRIAL 3
                        labTotalTauA[lab] += np.sum(labelImTau[...,lab]*weightsFine[...,None])
                    #labelImbin = labelImbin/labelImbinZ # normalize so all add up to 1 
                print("All should add up to 1 after first transition")
                print(np.unique(np.sum(labelImbin,axis=-1)))
                print("Labels before and after ")
                for lab in range(numLab):
                    print(labTotal[lab])
                    print(labTotalA[lab])
                print("Tau in Labels before and after ")
                for lab in range(numLab):
                    print(labTotalTau[lab])
                    print(labTotalTauA[lab])

                # save slice information as varifold (still in grid format)
                savenameBase = '/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/Block' + str(blockNum+1) + '/' + stainString + '/3DVarifold/VarifoldPartFiles/L' + str(locationNum) + '_FineDown' + str(exp) + '_jacobW' + str(weightBy) + '_' + str(brainNum) + 'only' + savename
                if (histoLab):
                    savenameBase = '/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/Block' + str(blockNum+1) + '/' + stainString + '/3DVarifold/VarifoldPartFiles/L' + str(locationNum) + '_FineDown' + str(exp) + '_2DLabels_jacobW' + str(weightBy) + '_' + str(brainNum) + 'only'
                feats = [histoTau,labelImTau,labelImbin,histoTauArea,histoTauOrient,histoTauRound]
                feats_dict = dict()
                for i in range(6):
                    feats_dict['feats'+str(i)] = feats[i]
                savenameTot = savenameBase + '.npz'
                np.savez(savenameTot,coords=histocoords_3D,weights=weightsFine,**feats_dict)
                convertToVarifoldExVivoPrelim(savenameBase)
                
                    # stack slice information onto others (keep as unraveled in terms of coordinates so each row denotes a particle)
                if (len(finePartCoords) < 1):
                    finePartCoords = np.reshape(histocoords_3D,(np.prod(histocoords_3D.shape[:-1]),histocoords_3D.shape[-1]))
                    finePartTau = np.ravel(histoTau)[...,None] # tau particles in tissue in fine particle (0 or 1)
                    finePartLab = np.reshape(labelImbin,(np.prod(labelImbin.shape[:-1]),labelImbin.shape[-1])) # vector valued fraction of area in region (0 or 1)
                    finePartTauLab = np.reshape(labelImTau,(np.prod(labelImbin.shape[:-1]),labelImbin.shape[-1])) # vector valued fraction of area in region (0 or 1)
                    finePartWeights = np.ravel(weightsFine)[...,None] # 0 or 0.00208 mm^2 (area of pixel)
                    finePartTauArea = np.ravel(histoTauArea)[...,None]
                    finePartTauOrient = np.ravel(histoTauOrient)[...,None]
                    finePartTauRound = np.ravel(histoTauRound)[...,None]
                else:
                    finePartCoords2 = np.reshape(histocoords_3D,(np.prod(histocoords_3D.shape[:-1]),histocoords_3D.shape[-1]))
                    finePartTau2 = np.ravel(histoTau)[...,None] # tau particles in tissue in fine particle (0 or 1)
                    finePartLab2 = np.reshape(labelImbin,(np.prod(labelImbin.shape[:-1]),labelImbin.shape[-1])) # vector valued fraction of area in region (0 or 1)
                    finePartTauLab2 = np.reshape(labelImTau,(np.prod(labelImbin.shape[:-1]),labelImbin.shape[-1])) # vector valued fraction of area in region (0 or 1)
                    finePartWeights2 = np.ravel(weightsFine)[...,None] # 0 or 0.00208 mm^2 (area of pixel)
                    finePartTauArea2 = np.ravel(histoTauArea)[...,None]
                    finePartTauOrient2 = np.ravel(histoTauOrient)[...,None]
                    finePartTauRound2 = np.ravel(histoTauRound)[...,None]

                    finePartCoords = np.vstack((finePartCoords,finePartCoords2))
                    finePartTau = np.vstack((finePartTau,finePartTau2)) # tau particles in tissue in fine particle (0 or 1)
                    finePartLab = np.vstack((finePartLab,finePartLab2)) # vector valued fraction of area in region (0 or 1)
                    finePartTauLab = np.vstack((finePartTauLab,finePartTauLab2)) # vector valued fraction of area in region (0 or 1)
                    finePartWeights = np.vstack((finePartWeights,finePartWeights2)) # 0 or 0.00208 mm^2 (area of pixel)
                    finePartTauArea = np.vstack((finePartTauArea,finePartTauArea2))
                    finePartTauOrient = np.vstack((finePartTauOrient,finePartTauOrient2))
                    finePartTauRound = np.vstack((finePartTauRound,finePartTauRound2))
    
    
    savenameTot = '/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/MultiBlock/' + stainString + '/3DVarifold/VarifoldPartFiles/All_FineDown' + str(exp) + '_jacobW' + str(weightBy) + '_' + str(brainNum) + 'only' + savename + '.npz'
    if (histoLab):
        savename = '/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/MultiBlock/' + stainString + '/3DVarifold/VarifoldPartFiles/All_FineDown' + str(exp) + '_2DLabels_jacobW' + str(weightBy) + '_' + str(brainNum) + 'only.npz'
    feats = [finePartTau,finePartTauLab,finePartLab,finePartTauArea,finePartTauOrient,finePartTauRound]
    feats_dict = dict()
    for i in range(6):
        feats_dict['feats'+str(i)] = feats[i]
    np.savez(savenameTot,coords=finePartCoords,weights=finePartWeights,**feats_dict)

    print("weights before and after")
    print(weightTotal)
    print(weightTotalA)
    
    print("Tau before and after (and sum)")
    print(tauTotal)
    print(tauTotalA)
    print(np.sum(labTotalTau))
    print(np.sum(labTotalTauA))
    
    print("Labels before and after ")
    for lab in range(numLab):
        print(labTotal[lab])
        print(labTotalA[lab])
    print("Tau in Labels before and after ")
    for lab in range(numLab):
        print(labTotalTau[lab])
        print(labTotalTauA[lab])
    
    return finePartCoords,feats,finePartWeights

def selectRegionOnly(filename,savenameBase,regions):
    '''
    Save new varifold where weights, total tau, tau in region/total weight, weight region / total weight, tau area, tau orient, tau round 
    are all with respect to total weights of region
    '''
    loc,weights,feats = loadVarifold(filename) # assume there are 6 features
    weightsNew = np.zeros_like(weights)
    
    tauLab = weights*feats[1][...,regions] # unnormalized
    tau = np.sum(tauLab,axis=-1)
    labFrac = weights*feats[2][...,regions] # unnormalized
    weightsNew = np.sum(labFrac,axis=-1)
    inds = weightsNew > 0
    
    recip = 1.0/weightsNew
    recip[recip == np.inf] = 0
    tauLab = tauLab*recip[...,None]
    tau = tau*recip
    labFrac = labFrac*recip[...,None]
    
    recipTau = np.squeeze(tau)/np.squeeze(feats[0])
    recipTau[recipTau == np.inf] = 1
    tauArea = np.squeeze(weights)*np.squeeze(feats[3])*recipTau*np.squeeze(recip) # assume uniform area so choose percent of tau that are taking
    tauOrient = feats[4] # just take as is
    tauRound = feats[5]
    print(tauOrient.shape)
    print(tauArea.shape)
    
        # only choose ones for which weightsNew > 0
    locNew = loc[inds]
    featsNew = [tau[inds],tauLab[inds],labFrac[inds],tauArea[inds],tauOrient[inds],tauRound[inds]]
    print(featsNew[0].shape)
    print(featsNew[1].shape)
    feats_dict = dict()
    for i in range(6):
        feats_dict['feats'+str(i)] = featsNew[i]
    np.savez(savenameBase+'.npz',coords=locNew,weights=weightsNew[inds],**feats_dict)
    convertToVarifoldExVivoPrelim(savenameBase)
    return

def selectRegionOnlyInter(filename,savenameBase,regions):
    '''
    Save new varifold where weights, total tau, tau in region/total weight, weight region / total weight, tau area, tau orient, tau round 
    are all with respect to total weights of region
    '''
    loc,weights,feats = loadVarifold(filename) # assume there are 6 features
    if (len(loc.shape)>2):
        loc = np.reshape(loc,(np.prod(loc.shape[0:-1]),3))
        weights = np.ravel(weights)
        weights = weights[...,None]
        for i in range(6):
            if (i == 1 or i == 2):
                feats[i] = np.reshape(feats[i],(np.prod(feats[i].shape[0:-1]),feats[i].shape[-1]))
            else:
                feats[i] = np.ravel(feats[i])[...,None]
                
    weightsNew = np.zeros_like(weights)
    
    tauLab = weights*feats[1][:,regions] # unnormalized
    tau = np.sum(tauLab,axis=-1)
    labFrac = weights*feats[2][:,regions] # unnormalized
    weightsNew = np.sum(labFrac,axis=-1)
    inds = weightsNew > 0
    
    recip = 1.0/weightsNew
    recip[recip == np.inf] = 0
    tauLab = tauLab*recip[...,None]
    tau = tau*recip
    labFrac = labFrac*recip[...,None]
    
    recipTau = np.squeeze(tau)/np.squeeze(feats[0])
    recipTau[recipTau == np.inf] = 1
    tauArea = np.squeeze(weights)*np.squeeze(feats[3])*recipTau*np.squeeze(recip) # assume uniform area so choose percent of tau that are taking
    tauOrient = feats[4] # just take as is
    tauRound = feats[5]
    print(tauOrient.shape)
    print(tauArea.shape)
    
        # only choose ones for which weightsNew > 0
    locNew = loc[inds]
    featsNew = [tau[inds],tauLab[inds],labFrac[inds],tauArea[inds],tauOrient[inds],tauRound[inds]]
    print(featsNew[0].shape)
    print(featsNew[1].shape)
    feats_dict = dict()
    for i in range(6):
        feats_dict['feats'+str(i)] = featsNew[i]
    np.savez(savenameBase+'.npz',coords=locNew,weights=weightsNew[inds],**feats_dict)
    convertToVarifoldExVivoPrelim(savenameBase)
    return

def selectZSection(filename,savenameBase,zLim):
    loc,weights,feats = loadVarifold(filename) # assume there are 6 features
    inds = (loc[...,2] > zLim[0])*(loc[...,2] < zLim[1]) # locations between region
    locNew = loc[inds]
    weightsNew = weights[inds]
    feats_dict = dict()
    for i in range(6):
        feats_dict['feats'+str(i)] = feats[i][inds]
        print(feats[i][inds].shape)
    np.savez(savenameBase+'.npz',coords=locNew,weights=weightsNew,**feats_dict)
    convertToVarifoldExVivoPrelim(savenameBase)
    return

def binToMulti(binLabIm):
    numLab = binLabIm.shape[-1]
    binLabRet = np.zeros_like(binLabIm[...,0])
    for l in range(numLab):
        binLabRet += binLabIm[...,0]*(l+1) # zero if is background
    return binLabRet

# Assume the image has all labels in it (get one for the background at end )
def multiToBin(multiLabIm, numLabs=0):
    if (numLabs == 0):
        numLabs = np.amax(multiLabIm) + 1 # include background
    print("max label in multi to bin is " + str(np.amax(multiLabIm)))
    print("using " + str(numLabs))
    binLabRet = np.zeros_like(multiLabIm)
    binLabRetList = []
    for i in range(numLabs):
        binLabRetList.append(binLabRet)
    # append one for background
    #binLabRetList.append(binLabRet)
    binLabRet = np.stack(binLabRetList,axis=-1)
    print("return label image shape")
    print(binLabRet.shape)
    for i in range(numLabs):
        binLabRet[...,i] = (multiLabIm == i+1).astype('float32')
    binLabRet[...,-1] = (multiLabIm == 0).astype('float32') # set last one to background
    print("should all add up to 1 in multiToBin (prints number that do not add up to 1):")
    print(np.sum(np.sum(binLabRet,axis=-1) != 1))
    return binLabRet

# Load Varifold
def loadVarifold(npzFile):
    fp = np.load(npzFile)
    fpLoc = fp['coords'].astype('float32')
    fpW = fp['weights'].astype('float32')
    fkeys = [i for i in fp.files if 'feats' in i]
    fpF = []
    for i in range(len(fkeys)):
        fpF.append(fp['feats'+str(i)].astype('float32'))
    return fpLoc,fpW,fpF

def loadVarifoldExVivo(npzFile):
    cp = np.load(npzFile)
    cpLoc0 = cp['X0s'].astype('float32')
    cpLoc1 = cp['X1s'].astype('float32')
    cpLoc2 = cp['X2s'].astype('float32')
    cpLoc = np.stack((cpLoc0,cpLoc1,cpLoc2),axis=-1)
    cpW = cp['cpWeight'].astype('float32')
    cpF = []
    cpF.append(cp['cpTau'].astype('float32')) # tau per area
    cpF.append(cp['cpLabTau'].astype('float32')) # tau in Lab per area
    cpF.append(cp['cpLab'].astype('float32')) # area in Lab per area
    cpF.append(cp['cpLabTauLab'].astype('float32')) # tau in Lab per area Lab
    return cpLoc,cpW,cpF

def convertToVarifoldExVivo(npzFileBase):
    cp = np.load(npzFileBase + '.npz')
    cpLoc = cp['coords'].astype('float32')
    X0s = cpLoc[...,0]
    X1s = cpLoc[...,1]
    X2s = cpLoc[...,2]
    cpWeight = cp['weights']
    fkeys = [i for i in cp.files if 'feats' in i]
    fpF = []
    for i in range(len(fkeys)):
        fpF.append(cp['feats'+str(i)].astype('float32'))
    cpTau = fpF[0]
    cpLabTau = fpF[1]
    cpLab = fpF[2]
    np.savez(npzFileBase+'_toVivo.npz',X0s=X0s,X1s=X1s,X2s=X2s,cpWeight=cpWeight,cpTau=cpTau,cpLabTau=cpLabTau,cpLab=cpLab)
    if (len(fpF) > 3):
        cpLabTauLab = fpF[3]
        np.savez(npzFileBase+'_toVivo.npz',X0s=X0s,X1s=X1s,X2s=X2s,cpWeight=cpWeight,cpTau=cpTau,cpLabTau=cpLabTau,cpLab=cpLab,cpLabTauLab=cpLabTauLab)
    # new files with std
    if (len(fpF) > 4):
        cpTauVar = fpF[4]
        cpLabTauLab2 = fpF[5]
        cpLabTauLab2Var = fpF[6]
        np.savez(npzFileBase+'_toVivo.npz',X0s=X0s,X1s=X1s,X2s=X2s,cpWeight=cpWeight,cpTau=cpTau,cpLabTau=cpLabTau,cpLab=cpLab,cpLabTauLab=cpLabTauLab,cpTauVar=cpTauVar,cpLabTauLab2=cpLabTauLab2,cpLabTauLab2Var=cpLabTauLab2Var)
    return

'''
Meant to be used for going to scaleGauss function in Varifold Framework
'''
def convertToVarifoldExVivoPrelim(npzFileBase):
    cp = np.load(npzFileBase + '.npz')
    cpLoc = cp['coords'].astype('float32')
    X0s = cpLoc[...,0]
    X1s = cpLoc[...,1]
    X2s = cpLoc[...,2]
    histocoords_3D = cpLoc
    weightsFine = cp['weights']
    print("weights shape")
    print(weightsFine.shape)
    fkeys = [i for i in cp.files if 'feats' in i]
    fpF = []
    for i in range(len(fkeys)):
        fpF.append(cp['feats'+str(i)].astype('float32'))
    histoTau = fpF[0]
    labelImTau = fpF[1]
    labelImbin = fpF[2]
    print("shapes of converted to ex vivo")
    print(histoTau.shape)
    print(labelImTau.shape)
    print(labelImbin.shape)
    np.savez(npzFileBase+'_toVivo.npz',histocoords_3D=histocoords_3D,weightsFine=weightsFine,histoTau=histoTau,labelImTau=labelImTau,labelImbin=labelImbin)
    return


def saveAsMat(weights,feats,savename,maxArea=None):
    dictToSave = dict()
    dictToSave['weights'] = weights
    if (len(weights.shape) < 2):
        dictToSave['weights'] = weights[...,None]
    if (maxArea is not None):
        dictToSave['maxArea'] = maxArea
    for i in range(len(feats)):
        if (len(feats[i].shape)<2):
            dictToSave['feats'+str(i)] = feats[i][...,None]
        else:
            dictToSave['feats'+str(i)] = feats[i]
    
    sp.io.savemat(savename,dictToSave)
    return

################################################################
# Fine to Coarse

# Main Module for Downsampling Varifold Data
def fineToCoarse(coarsePartLoc,finePartLoc,finePartFeat,finePartWeight,sKern,sKernParams, fKern,fKernParams,conserveMem=False):
    '''
    Args:
    coarsePartLoc = *x3 matrix with x,y,z locations of coarse particles (i.e. Nx3 is a list of irregularly sampled)
    finePartLoc = *x3 matrix with x,y,z locations of fine particles
    finePartFeat = [*xF] list of matrices with feature values per feature
    finePartWeight = * matrix with fine particle weights
    sKern = spatial kernel to apply to adjusted weights 
    fKern = feature kernel that takes in all dimensions of features 
    sKernParams = list of parameters to use with sKern
    fKernParams = list of parameters to use with fKern
    cpdimF = dimension of feature space out (for preallocation of memory)
    
    Returns:
    cpWeights = Weights of coarse particles in same format as locations
    cpFeatures = vector of feature values for each coarse particle
    '''
    
    cpWeights = np.zeros(coarsePartLoc.shape[:-1])
    cpFeatures = [] # list of features
    finePartFeatRav = [] # list of unraveled features (fp x fv)
    numFP = np.prod(finePartLoc.shape[:-1])
    numCP = np.prod(coarsePartLoc.shape[:-1])
    print("number of fine and coarse particles: ")
    print(numFP)
    print(numCP)
    for f in finePartFeat:
        if (np.prod(f.shape) == numFP):
            finePartFeatRav.append(f.ravel())
        else:
            finePartFeatRav.append(np.reshape(f,(np.prod(f.shape[:-1]),f.shape[-1]))) # unraveled feature values (fp x fv) per feature
    numFPfeat = len(finePartFeat)
    print("number of features in rav: " + str(len(finePartFeatRav)))
    
    # assume kernel function is a conserve memory one that returns file paths of fine particle info
    if (conserveMem):
        cpWeightsRav, newFineWeightsRav = sKern(coarsePartLoc,finePartLoc,finePartWeight,sKernParams) #list of unraveled fine particle weights per unraveled coarse particles (cp x fp)
        print('finished spatial kernel')
        #cpWeightsRav = sum(np.sum(k,axis=-1) for k in newFineWeightsRav)
        '''
        cpWeightsRav = np.zeros(numCP)
        for nFN in newFineWeightsRav:
            fpInfo = np.load(nFN)['wA'].astype('float32')
            print(fpInfo.shape)
            cpWeightsRav += np.sum(fpInfo,axis=-1)
        '''
        cpWeights = np.reshape(cpWeightsRav,cpWeights.shape) # check to make sure this reverses correctly
        reWeights = 1./cpWeights
        reWeights[reWeights==np.inf] = 0 # 0 or 1/w_j if w_j \neq \infty
        reWeightsRav = 1./cpWeightsRav
        reWeightsRav[reWeightsRav == np.inf] = 0
        print("After space kernel")
        print(np.amax(cpWeightsRav))
        print(np.amax(newFineWeightsRav))
        
        cpFeaturesRav = []
        for j in range(numCP): # size fp x 0
            if (j % 10 == 0):
                print("On CP " + str(j))
            #cpw = reWeightsRav[j]*np.hstack((k[j,:].toarray() for k in newFineWeightsRav)) # assuming sparse representation
            
            cpw = np.asarray([])
            firstOne = True
            for nFN in newFineWeightsRav:
                fpInfo = np.load(nFN)['wA'].astype('float32')
                if (firstOne):
                    cpw = fpInfo[j,:]
                    firstOne = False
                else:
                    cpw = np.hstack((cpw,fpInfo[j,:]))
            
            cpw = reWeightsRav[j]*cpw
            print("cpw shape is " + str(cpw.shape))
            
            nu_j = [] # feature varifolds as intermediate (sum of particles), list of 1xfpxfv
            for f in finePartFeatRav:
                if (np.prod(f.shape) == numFP):
                    a = cpw*f
                    nu_j.append(a)
                else:
                    a = cpw[...,None]*f
                    nu_j.append(a) # separate matrix per feature (1 x fp x fv) --> cp x fp x fv 

            cpF_j = fKern(nu_j,fKernParams) # returns list of arrays of fv for single CP
            if (j == 0):
                # add np.zeros equal to #cp x fv for each features
                numCPFeat = len(cpF_j)
                for cpF in cpF_j:
                    cpFeaturesRav.append(np.zeros((numCP,np.prod(cpF.shape)))) # assume all 1s or 0s other than fv
            for k in range(numCPFeat):
                cpFeaturesRav[k][j] = cpF_j[k]
            j += 1

        for f in cpFeaturesRav:
            print("Feature values are")
            print(f)
            if (np.prod(f.shape) == numCP):
                cpFeatures.append(np.reshape(f,cpWeights.shape))
            else:
                cpFeatures.append(np.reshape(f,cpWeights.shape + (f.shape[-1],)))
        
    else:
        newFineWeightsRav = sKern(coarsePartLoc,finePartLoc,finePartWeight,sKernParams) #unraveled fine particle weights per unraveled coarse particles (cp x fp)
        print("number of fine particles is " + str(newFineWeightsRav.shape[1]))
        cpWeightsRav = np.sum(newFineWeightsRav,axis=-1)
        cpWeights = np.reshape(cpWeightsRav,cpWeights.shape) # check to make sure this reverses correctly
        print("cp weights are ")
        print(cpWeights)
        reWeights = 1.0/cpWeights
        reWeights[reWeights==np.inf] = 0.0 # 0 or 1/w_j if w_j \neq \infty
        reWeightsRav = 1.0/cpWeightsRav
        reWeightsRav[reWeightsRav==np.inf] = 0.0
        print(reWeightsRav)

        '''
        For computation reasons, I make the feature kernels only take in 1 cp info at a time
        '''

        w = newFineWeightsRav*reWeightsRav[...,None] # cp x fp * cp x 1 = cp x fp
        j = 0
        cpFeaturesRav = []
        for cpw in w: # size fp x 0, equal to w_i k(x_i, y_j) / w_j (j is row, all i in each cpw)
            nu_j = [] # feature varifolds as intermediate (sum of particles), list of 1xfpxfv
            for f in finePartFeatRav:
                if (np.prod(f.shape) == numFP):
                    a = cpw*f 
                    a = f # for new fKern
                    nu_j.append(a)
                else:
                    a = cpw[...,None]*f
                    a = f # for new fKern
                    nu_j.append(a) # separate matrix per feature (1 x fp x fv) --> cp x fp x fv 

            #cpF_j = fKern(nu_j,fKernParams) # returns list of arrays of fv for single CP
            cpF_j = fKern(cpw,nu_j,fKernParams)
            if (j == 0):
                # add np.zeros equal to #cp x fv for each features
                numCPFeat = len(cpF_j)
                for cpF in cpF_j:
                    cpFeaturesRav.append(np.zeros((numCP,np.prod(cpF.shape)))) # assume all 1s or 0s other than fv
            for k in range(numCPFeat):
                cpFeaturesRav[k][j] = cpF_j[k]
            j += 1

        for f in cpFeaturesRav:
            print("Feature values are")
            print(f)
            if (np.prod(f.shape) == numCP):
                cpFeatures.append(np.reshape(f,cpWeights.shape))
            else:
                cpFeatures.append(np.reshape(f,cpWeights.shape + (f.shape[-1],)))

    return cpWeights,cpFeatures


def scaleGauss(coarseGridNPZ,brainNum,blocks,sigma,numLab,stainNum,numDown,intdenom,gauss,byLabel,subsetSuff,regions='',toVivo=False):
    '''
    Args:
        1) coarseGridNPZ = output locations of coarse particles
        2) brainNum, blocks
        3) sigma = for gaussian (and used in coarse grid sampling where sampling is assumed) sigma/intdenom
        4) numLab = assumed number of labels in label function (INCLUDES BACKGROUND = UNLABELED TISSUE)
        5) numDown = number of times fine grid was downsampled in first step / translation to varifold (no downsampling = 0)
        6) gauss = True denotes use of gaussian rather than Haar
        7) byLabel = True denotes storing additional values by label 

    Notes:
        Feature kernel is taken to be identity for all features here (i.e. assumes feature values in desired form as input)
        Weights = area sampled, Tau = particles / mm^2, Labels = fraction of area sampled in designated region
    '''
    if (stainNum == 0):
        stainString = "Amyloid"
    elif (stainNum == 2):
        stainString = "Tau"
    
    if (subsetSuff is None):
        subsetSuff = ''
        
    cgParams = np.load(coarseGridNPZ)
    x0samp = cgParams['x0samp'].astype('float32')
    x1samp = cgParams['x1samp'].astype('float32')
    x2samp = cgParams['x2samp'].astype('float32')
    X0s = cgParams['X0s'].astype('float32')
    X1s = cgParams['X1s'].astype('float32')
    X2s = cgParams['X2s'].astype('float32')
    print(x0samp)
    print(x1samp)
    print(x2samp)
    
    if (len(sigma) > 1):
        sigmaString = str(sigma[0]) + str(sigma[1]) + str(sigma[2])
        sigma0 = sigma[0]
        sigma1 = sigma[1]
        sigma2 = sigma[2]
    else:
        sigmaString = str(sigma[0])
        sigma0 = sigma[0]
        sigma1 = sigma[0]
        sigma2 = sigma[0]
    
    # Course Particles
    cpWeight = np.zeros_like(X0s).astype('float32') # will be transparency (sum of weights of fine particles)
    cpTau = np.zeros_like(X0s).astype('float32') # will be expected value of tau particle based on empirical distribution of fine particles
    cpLab = np.zeros((X0s.shape[0],X0s.shape[1],X0s.shape[2],numLab)).astype('float32') # expected value of each indicator based on empirical distribution of fine particles; background = last channel
    cpLabTau = None
    cpLabTauLab = None
    if (byLabel):
        cpLabTau = np.zeros_like(cpLab)
        cpLabTauLab = np.zeros_like(cpLab)
    
    weightTotal = 0 # sum weights of fine particles
    tauTotal = 0 # total particles (sum of fine particle tau function * fine particle initial weight)
    labTotal = np.zeros((numLab,1))
    labTotalTau = np.zeros((numLab,1))
    epsilon = 2e-10
    #epsilon = 0.0
    
    for blockNum in blocks:
        if (brainNum == 2):
            high = 15
            start = 0
            slices = [5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]
            if (blockNum == 2):
                high = 7
        elif (brainNum == 5):
            slices = [1,5,10,15,20,25,30,35,40,45,50,55,60,65,70]
            start = 0
            if (blockNum == 0):
                high = 13
            elif (blockNum == 1):
                high = 14
            elif (blockNum == 2):
                slices = [3,7,12,17,22,27,32,37]
                high = 8
        for sl in range(start,high):
            locationNum = slices[sl]
            print("On block " + str(blockNum+1) + " and loc " + str(locationNum))
            savename = '/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/MultiBlock' + '/' + stainString + '/3DVarifold/Slices/3D_intermediate_down' + str(numDown) + '_block' + str(blockNum+1) + '_slice' + str(sl) + '_newmai' + subsetSuff + '.npz'
            if (toVivo):
                savename = '/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/Block' + str(blockNum+1) + '/' + stainString + '/3DVarifold/VarifoldPartFiles/L' + str(locationNum) + regions + '_FineDown' + str(numDown) + subsetSuff + '_toVivo.npz' # subsetSuff for toVivoshould be W2 or W3
            if (not path.exists(savename)):
                continue
            params = np.load(savename)
            histocoords_3D = params['histocoords_3D'].astype('float32')
            histoTau = params['histoTau'].astype('float32')
            weightsFine = params['weightsFine'].astype('float32')
            labelImbin = params['labelImbin'].astype('float32') # fraction of area in region
            labelImTau = params['labelImTau'].astype('float32') # number of Tau particles/labelImbin*weightsFine in region
            
            weightTotal += np.sum(weightsFine)
            tauTotal += np.sum(weightsFine * histoTau)
            for l in range(numLab):
                labTotal[l] += np.sum(np.squeeze(labelImbin[...,l])*weightsFine) # sum of the area of given label
                labTotalTau[l] += np.sum(labelImTau[...,l]*weightsFine[...,None])
            histoTauRav = np.ravel(histoTau[...,None])
            histoTauRav = histoTauRav[...,None].T
            # Calculate weights for fine particles in this slice only
            #totNumParts = histocoords_3D.shape[0]*histocoords_3D.shape[1]*histocoords_3D.shape[2] # last dimension should be 1
            totNumParts = np.prod(histocoords_3D.shape[0:-1]) # can have raveled 
            at1 = np.ravel(histocoords_3D[...,0])
            at2 = np.ravel(histocoords_3D[...,1])
            at3 = np.ravel(histocoords_3D[...,2])
            at4 = np.reshape(labelImbin,(np.prod(labelImbin.shape[0:-1]),labelImbin.shape[-1]))
            if (byLabel):
                at5 = np.reshape(labelImTau,(np.prod(labelImTau.shape[0:-1]),labelImTau.shape[-1]))
                at5 = at5.T
            at1 = at1[...,None].T        
            at2 = at2[...,None].T
            at3 = at3[...,None].T
            at4 = at4.T
            weightsRav = np.ravel(weightsFine[...,None])
            weightsRav = weightsRav[...,None].T

            # Kernel for calculating weights (summing over J)
            if (gauss):
                zeta = np.zeros_like(at1).astype('float32') # Normaliization so all info is included 
               
                cpLabSlice = np.zeros((X0s.shape[0],X0s.shape[1],X0s.shape[2],numLab)).astype('float32') # independent of label now
                cpTauSlice = np.zeros_like(X0s).astype('float32')
                if (byLabel):
                    cpLabTauSlice = np.zeros_like(cpLabSlice).astype('float32')
                    cpLabTauLabSlice = np.zeros_like(cpLabSlice).astype('float32')
                else:
                    cpLabTauSlice = None
                    cpLabTauLabSlice = None
                
                '''
                Loop for calculating zeta alone
                '''
                pointsPos = []
                for x0 in range(x0samp.shape[0]):
                    x0a = x0samp[x0]
                    summ0 = np.exp((1.0/(2.0*sigma0**2))*(-(x0a - at1)**2))
                    if (np.prod(at1.shape) < 1 or np.max(summ0) < epsilon):
                        continue # no datapoints close enough here
                    for x1 in range(x1samp.shape[0]):
                        x1a = x1samp[x1]
                        summ01 = summ0*np.exp((1.0/(2.0*sigma1**2))*(-(x1a - at2)**2))
                        if (np.max(summ01) < epsilon):
                            continue
                        for x2 in range(x2samp.shape[0]):
                            x2a = x2samp[x2]
                            #print("on " + str(x0) + ", " + str(x1) + ", " + str(x2))
                            #summ = (1.0/(np.sqrt(2.0*np.pi*sigma**2))**3)*np.exp((1.0/(2.0*sigma**2))*(-(x0a - at1)**2 - (x1a - at2)**2 - (x2a - at3)**2)) # zero out ones not in label
                            summ = summ01*(1.0/(np.sqrt((2.0*np.pi)**3*sigma0**2*sigma1**2*sigma2**2)))*np.exp((1.0/(2.0*sigma2**2))*(-(x2a - at3)**2))
                            zeta += summ
                            pointsPos.append((x0,x1,x2))

                
                kde3D = np.zeros_like(X0s).astype('float32')

                zetaD = 1.0/zeta
                zetaD[zetaD==np.inf] = 0.0
                print("number of particles with nonzero weight that had 0 contribution: " )
                print(np.sum((zetaD == 0.0)*(weightsRav > 0)))
                
                for p in pointsPos:
                    x0a = x0samp[p[0]]
                    x1a = x1samp[p[1]]
                    x2a = x2samp[p[2]]
                    summ0 = np.exp((1.0/(2.0*sigma0**2))*(-(x0a - at1)**2))
                    summ01 = summ0*np.exp((1.0/(2.0*sigma1**2))*(-(x1a - at2)**2))
                    summ = summ01*(1.0/(np.sqrt((2.0*np.pi)**3*sigma0**2*sigma1**2*sigma2**2)))*np.exp((1.0/(2.0*sigma2**2))*(-(x2a - at3)**2))
                    summ = summ*weightsRav
                    summ = summ*zetaD #now gives precise weights
                    kde3D[p[0],p[1],p[2]] = np.sum(summ,axis=-1)
                    cpTauSlice[p[0],p[1],p[2]] = np.sum((summ*histoTauRav),axis=-1)
                    for lab in range(numLab):
                        cpLabSlice[p[0],p[1],p[2],lab] = np.sum((summ*at4[lab,...]),axis=-1)
                        if (byLabel):
                            cpLabTauSlice[p[0],p[1],p[2],lab] = np.sum((summ*at5[lab,...]),axis=-1)
                    
                cpWeight += kde3D
                cpTau += cpTauSlice
                cpLab += cpLabSlice
                cpLabTau += cpLabTauSlice

                # Save Intermediate slice only information
                savename3D = '/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/MultiBlock' + '/' + stainString + '/3DVarifold/Slices/3D_coarseInfo_sigma' + str(sigmaString) + '_intdenom' + str(intdenom) + '_gauss' + str(gauss) + '_block' + str(blockNum+1) + '_slice' + str(sl) + subsetSuff + '.npz'
                np.savez(savename3D,kde3D=kde3D,tauSlice=cpTauSlice,labSlice=cpLabSlice,labTauSlice=cpLabTauSlice)
            
            # Kernel of voxel only 
            else:
                kdeSlice = np.zeros_like(X0s) # independent of label now
                cpLabSlice = np.zeros((X0s.shape[0],X0s.shape[1],X0s.shape[2],numLab)) # independent of label now
                cpTauSlice = np.zeros_like(X0s)


                for x0 in range(x0samp.shape[0]):
                    x0a = x0samp[x0]
                    summ0 = (x0a - 0.5*(sigma/intdenom) <= at1)*(at1 < x0a + 0.5*(sigma/intdenom))
                    if (np.sum(summ0 != 0) < 1):
                        continue
                    for x1 in range(x1samp.shape[0]):
                        x1a = x1samp[x1]
                        summ01 = summ0*(x1a - 0.5*(sigma/intdenom) <= at2)*(at2 < x1a + 0.5*(sigma/intdenom))
                        if (np.sum(summ01 != 0) < 1):
                            continue
                        for x2 in range(x2samp.shape[0]):
                            x2a = x2samp[x2]
                            print("on " + str(x0) + ", " + str(x1) + ", " + str(x2))
                            #summ = (1.0/(np.sqrt(2.0*np.pi*sigma**2))**3)*np.exp((1.0/(2.0*sigma**2))*(-(x0a - at1)**2 - (x1a - at2)**2 - (x2a - at3)**2)) # zero out ones not in label
                            summ = summ01*(x2a - 0.5*(sigma/intdenom) <= at3)*(at3 < x2a + 0.5*(sigma/intdenom)) 
                            summ = summ*weightsRav
                            #summ = summ*weightsRav # should give exact weights per fine particles
                            cpWeight[x0,x1,x2] += np.sum(summ) # sum of all area
                            kdeSlice[x0,x1,x2] = np.sum(summ)
                            for lab in range(numLab):
                                summLab = (at4[lab,...])*summ # won't be a binary value, unnormalized weights
                                cpLab[x0,x1,x2,lab] += np.sum(summLab)
                                cpLabSlice[x0,x1,x2,lab] = np.sum(summLab)
                            summTau = (histoTauRav * summ)
                            cpTau[x0,x1,x2] += np.sum(summTau)
                            cpTauSlice[x0,x1,x2] = np.sum(summTau)
                
                # Save Intermediate slice only information
                savename3D = '/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/MultiBlock' + '/' + stainString + '/3DVarifold/Slices/3D_coarseInfo_sigma' + str(sigmaString) + '_intdenom' + str(intdenom) + '_gauss' + str(gauss) + '_block' + str(blockNum+1) + '_slice' + str(sl) + regions + subsetSuff + '.npz'
                np.savez(savename3D,kde3D=kdeSlice,tauSlice=cpTauSlice,labSlice=cpLabSlice)
            

    # Normalize everything
    # No other normalization because here, assumed that quantities already put in particles / mm^2 
    normWeights = 1./cpWeight
    normWeights[normWeights==np.inf] = 0
    cpTau = normWeights*cpTau
    print("Labels before and after")
    for lab in range(numLab):
        print(labTotal[lab])
        print(np.sum(cpLab[...,lab]))
        cpLab[...,lab] = cpLab[...,lab]*normWeights
    if (byLabel):
        print("Tau per Label before and after")
        totalSumOfTau = 0
        for lab in range(numLab):
            print(labTotalTau[lab]) # should be number of tau particles per region
            cpLabTau[...,lab] = cpLabTau[...,lab]*normWeights # TRIAL 3 -- keep as identity
            print(np.sum(cpLabTau[...,lab]*cpWeight))
            totalSumOfTau += np.sum(cpLabTau[...,lab]*cpWeight)
            recipLab = 1.0/cpLab[...,lab]
            recipLab[recipLab == np.inf] = 0
            cpLabTauLab[...,lab] = cpLabTau[...,lab]*recipLab # tau per mm of region 
    print("weights before and after")
    print(weightTotal)
    print(np.sum(cpWeight))
    
    print("Tau before and after and sum")
    print(tauTotal)
    print(np.sum(cpWeight*cpTau))
    print(totalSumOfTau)

    # Save variables to debug
    savename3D = '/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/MultiBlock' + '/' + stainString + '/3DVarifold/3D_coarseInfo_sigma' + str(sigmaString) + '_intdenom' + str(intdenom) + '_gauss' + str(gauss) + regions + subsetSuff + '.npz'
    np.savez(savename3D,cpWeight=cpWeight,cpTau=cpTau,cpLab=cpLab,cpLabTau=cpLabTau,cpLabTauLab=cpLabTauLab,X0s=X0s,X1s=X1s,X2s=X2s)
    return

def getCoarseGrid(brainNum,blocks,stainNum,sigma,intdenom,dims=None,minX=None,maxX=None,saveName=None,subsetSuff=''):
    '''
    dims = list of directions (dimensions) along which to sample points if reducing to 2D or 3D.
    If particular direction is not included in list, use average of min and max 
    '''
    
    if (stainNum == 0):
        stainString = "Amyloid"
    elif (stainNum == 2):
        stainString = "Tau"

    if (len(sigma) > 1):
        sigmaString = str(sigma[0]) + str(sigma[1]) + str(sigma[2])
        sigma0 = sigma[0]
        sigma1 = sigma[1]
        sigma2 = sigma[2]
    else:
        sigmaString = str(sigma[0])
        sigma0 = sigma[0]
        sigma1 = sigma[0]
        sigma2 = sigma[0]
    
    if (minX is None):
        minX0 = 10000000
        maxX0 = -100000000
        minX1 = 100000000
        maxX1 = -100000000
        minX2 = 100000000
        maxX2 = -100000000


        # find min and max out of ALL blocks
        for blockNum in blocks:
            high = 15
            if (brainNum == 2 and blockNum == 2):
                high = 7
            for sl in range(0,high):
                locationNum = sl+5
                histocoords_3D_npz = np.load('/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/2DHisto_to_MRI_coords_block' + str(blockNum+1) + '_stain' + str(stainNum) + '_slice' + str(sl) + subsetSuff + '.npz')
                histocoords_3D = histocoords_3D_npz['coordsT_mai']
                histocoords_3D = histocoords_3D.astype('float32')
                histocoords_3D_npz.close()

                if (np.amin(histocoords_3D[...,0]) < minX0):
                    minX0 = np.amin(histocoords_3D[...,0])
                if (np.amax(histocoords_3D[...,0]) > maxX0):
                    maxX0 = np.amax(histocoords_3D[...,0])
                if (np.amin(histocoords_3D[...,1]) < minX1):
                    minX1 = np.amin(histocoords_3D[...,1])
                if (np.amax(histocoords_3D[...,1]) > maxX1):
                    maxX1 = np.amax(histocoords_3D[...,1])
                if (np.amin(histocoords_3D[...,2]) < minX2):
                    minX2 = np.amin(histocoords_3D[...,2])
                if (np.amax(histocoords_3D[...,2]) > maxX2):
                    maxX2 = np.amax(histocoords_3D[...,2])
            print("Finished block " + str(blockNum))

        # Set min and max to be 1 less or 1 more than min / max rounded to integer
        minX0 = np.floor(minX0) - 1
        maxX0 = np.ceil(maxX0) + 1
        minX1 = np.floor(minX1) - 1
        maxX1 = np.ceil(maxX1) + 1
        minX2 = np.floor(minX2) - 1
        maxX2 = np.ceil(maxX2) + 1
    else:
        minX0 = minX[0]
        minX1 = minX[1]
        minX2 = minX[2]
        maxX0 = maxX[0]
        maxX1 = maxX[1]
        maxX2 = maxX[2]
        
    if (dims is None):
        x0samp = np.arange(minX0,maxX0,sigma0/intdenom)
        x1samp = np.arange(minX1,maxX1,sigma1/intdenom)
        x2samp = np.arange(minX2,maxX2,sigma2/intdenom)
    else:
        if 0 in dims:
            x0samp = np.arange(minX0,maxX0,sigma0/intdenom)
        else:
            x0samp = np.asarray([(maxX0+minX0)/2.0])
        if 1 in dims:
            x1samp = np.arange(minX1,maxX1,sigma1/intdenom)
        else:
            x1samp = np.asarray([(maxX1+minX1)/2.0])
        if 2 in dims:
            x2samp = np.arange(minX2,maxX2,sigma2/intdenom)
        else:
            x2samp = np.asarray([(maxX2+minX2)/2.0])
            # use below if want multi
            #x2samp = np.asarray([(maxX2+minX2)/4.0,(maxX2+minX2)/2.0,3.0*(maxX2+minX2)/4.0])
           
                
    X0s,X1s,X2s = np.meshgrid(x0samp,x1samp,x2samp,indexing='ij') # coarse particle locations
    if (saveName is None):
        savename = '/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/MultiBlock' + '/' + stainString + '/3DVarifold/3D_coarseGrid_sigma' + str(sigmaString) + '_int' + str(intdenom)+'.npz'
    else:
        #savename = '/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/MultiBlock' + '/' + stainString + '/3DVarifold/3D_coarseGrid_sigma' + str(sigmaString) + '_int' + str(intdenom)+'_newmai' + subsetSuff + '.npz'
        savename = saveName
    np.savez(savename, x0samp=x0samp,x1samp=x1samp,x2samp=x2samp,X0s=X0s,X1s=X1s,X2s=X2s)
    
    return





################################################################
# Kernels for saving Fine to Coarse 

def Id(nu_jW,nu_jF,params=None):
    # sum up each of the particles together 
    '''
    mu = []
    for f in nu_j:
        if (np.prod(f.shape) == f.shape[0]):
            el = np.zeros((1,1))
            el[0,0] = np.sum(f,axis=0)
        else:
            el = np.sum(f,axis=0)
        mu.append(el) # collapse the fine particles = fpxfv
    return mu # assume nu in cp x fp x f for each feature
    '''

    mu = []
    for n in range(len(nu_jF)): # each coarse particle
        if (np.prod(nu_jF[n].shape) == nu_jF[n].shape[0]):
            f = nu_jF[n]*nu_jW
            el = np.zeros((1,1))
            el[0,0] = np.sum(f,axis=0)
        else:
            f = nu_jF[n]*nu_jW[...,None]
            el = np.sum(f,axis=0)
        mu.append(el)
    return mu



def nearestNeighbor(coarsePartLoc,finePartLoc,finePartWeight,sparams):
    '''
    Assigns Points to nearest (weights will be 0 or 1)  
    
    Args:
    coarsePartLoc = *x3 matrix with x,y,z locations of coarse particles
    finePartLoc = *x3 matrix with x,y,z locations of fine particles
    finePartWeight = * matrix with fine particle weights
    
    sparams = savelists as
    Returns fine particle adjusted weights per coarse particle in numcpxnumpfinepart matrix
    '''
    savenpy = sparams[0]
    finePW = np.ravel(finePartWeight)
    print("fine particle weight is")
    print(finePartWeight.shape)
    at1 = np.ravel(finePartLoc[...,0]).astype('float32')
    at2 = np.ravel(finePartLoc[...,1]).astype('float32')
    at3 = np.ravel(finePartLoc[...,2]).astype('float32')
    
    weightsRav = np.ravel(finePartWeight[...,None]).astype('float32')
    weightsRav = weightsRav[...,None].T
    
    x0samp = np.ravel(coarsePartLoc[...,0]).astype('float32') # cp
    x1samp = np.ravel(coarsePartLoc[...,1]).astype('float32') # cp
    x2samp = np.ravel(coarsePartLoc[...,2]).astype('float32') # cp
    numCP = x0samp.shape[0]
    numFP = at1.shape[0]
    print(numCP)
    print(numFP)
    
    at = np.vstack((at1,at2,at3)) #3xfp
    print("should be 3xfp")
    print(at.shape)
    cpAssoc = np.zeros((1,numFP)) - 1
    minDist = np.zeros((1,numFP)) + 1000000.0

    pointsPos = [] # coarse particles with nonzero contribution from fine particles in resampling
    for iCP in range(numCP):
        x0a = x0samp[iCP]
        x1a = x1samp[iCP]
        x2a = x2samp[iCP]
        
        xa = np.asarray([x0a,x1a,x2a])
        xa = xa[...,None]

        diff = xa-at # gives # fine parts by number of fine parts but only need diag
        diff2 = diff[0,...]**2 + diff[1,...]**2 + diff[2,...]**2 # only if all three are below do we capture it
        cpAssoc[diff2 < minDist] = iCP
        minDist = diff2*(diff2<minDist) + minDist*(diff2>=minDist)
    weightsAdj = np.zeros((numCP,numFP)).astype('float32')
    minDistTot = np.zeros((numCP,numFP)).astype('float32')
    for iCP in range(numCP):
        weightsAdj[iCP,:] = (cpAssoc == iCP).astype('float32')*finePW[None,...]
        minDistTot[iCP,:] = (cpAssoc == iCP).astype('float32')*minDist
    print("should be all accounted for")
    print(np.sum(cpAssoc < 0))
    np.save(savenpy,minDistTot)
            
    return weightsAdj

def interToCoarseExVivoHuman11TnonLinear(nu_jW,nu_jF,params=None):
    '''
    Similar to interToCoarseExVivoHuman11T except nu_j is full varifold representation (weights and features separate)
    nu_jF = list of features (one line per particle)
    nu_jW = list of NORMALIZED fine particle weights (i.e. w_i * k(x_i,y_j)/w_j)
    
    From downsampleExVivoAndConvert:
    
    '''
    
    mu = []
    for n in range(len(nu_jF)): # each coarse particle
        if ((len(nu_jF) == 7 or len(nu_jF) == 4) and n == 3):
            mu.append(np.zeros((1,1)))
            continue
        elif (n==3):
            mu.append(np.zeros((1,1))) # add in 

        if (np.prod(nu_jF[n].shape) == nu_jF[n].shape[0]):
            f = nu_jF[n]*nu_jW
            el = np.zeros((1,1))
            el[0,0] = np.sum(f,axis=0)
        else:
            f = nu_jF[n]*nu_jW[...,None]
            el = np.sum(f,axis=0)
        mu.append(el)
    if (nu_jF[1].shape != nu_jF[2].shape):
        print("shapes are not equal as expected")
    
    # add tau in region / region mm^2
    # Based on Old Way (i.e. E[tau^2]/E[labelregion^2])
    if (len(nu_jF[2].shape) < 2 or nu_jF[2].shape[-1] == 1):
        extraD = 1.0/np.sum(nu_jF[2][...,None]*nu_jW[...,None],axis=0)
        extraD[extraD==np.inf]=0
        extra = np.sum(nu_jF[1][...,None]*nu_jW[...,None],axis=0)*extraD
    else:
        extraD = 1.0/np.sum(nu_jF[2]*nu_jW[...,None],axis=0)
        extraD[extraD==np.inf]=0
        extra = np.sum(nu_jF[1]*nu_jW[...,None],axis=0)*extraD
    #mu.append(extra)
    mu[3] = extra
    
    # add std tau / mm^2 as E[(tau/mm^2)^2] - E[tau/mm^2]^2 (both expected values taken with respect to empirical)
    var = nu_jW*nu_jF[0]*nu_jF[0] # E[squared value]
    el = np.zeros((1,1))
    el[0,0] = np.sum(var,axis=0)
    el[0,0] = el[0,0] - (mu[0]*mu[0])
    mu.append(el)
    
    wMult = nu_jW[...,None]
    if (np.prod(nu_jF[2].shape) == np.prod(nu_jW.shape)):
        extraD = 1.0/nu_jF[2][...,None] 
        extraD[extraD==np.inf]=0
        extra = np.sum(wMult*nu_jF[1][...,None]*extraD,axis=0)
        mu.append(extra)
    
        # E[(tau/region)^2] - (E[tau/region])^2
        extraSq = np.sum(extraD*extraD*nu_jF[1][...,None]*nu_jF[1][...,None]*wMult,axis=0)
        extraSq = extraSq - (extra*extra)
        mu.append(extraSq)
    else:
        # E[tau/region]
        extraD = 1.0/nu_jF[2] 
        extraD[extraD==np.inf]=0
        extra = np.sum(wMult*nu_jF[1]*extraD,axis=0)
        mu.append(extra)

        # E[(tau/region)^2] - (E[tau/region])^2
        extraSq = np.sum(extraD*extraD*nu_jF[1]*nu_jF[1]*wMult,axis=0)
        extraSq = extraSq - (extra*extra)
        mu.append(extraSq)
    return mu


###############################################################
# Wrapper Functions

def selectMTL(brainNum,stain,mriSuffix,dateNew):
    '''
    Select regions: 
        - amygdala
        - hata
        - amygdala + hata
        - erc + extension
        - ca1
        - subiculum
        - hippocampus All (ca1, ca2, ca3, subiculum, presubiculum, parasubiculum, hilus, molecular, granular)
    '''
    savename = '_'+mriSuffix+'_'+dateNew
    filename = '/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/MultiBlock/' + stain + '/3DVarifold/VarifoldPartFiles/All_FineDown4_jacobW2_' + str(brainNum) + 'only' + savename + '.npz'
    if (brainNum == 2):
        amy = [2]
        hata = [0]
        amyHata = [0,2]
        ercExt = [6,7]
        ca1 = [3]
        sub = [13]
        hippo = [3,4,5,8,9,10,11,12,13]
        hippoSub = [3,4,5,8,9,10,13]
        hippoAlv = [3,4,5,8,9,10,11,12,13,1]
        hippoSubAlv = [3,4,5,8,9,10,13,1]
    elif (brainNum == 5):
        amy = [14]
        hata = [7]
        amyHata = [7,14]
        ercExt = [5,13]
        ca1 = [1]
        sub = [12]
        hippo = [1,2,3,4,6,8,9,10,11,12]
        hippoSub = [1,2,3,4,6,8,9,12]
        hippoAlv = [1,2,3,4,6,8,9,10,11,12,0]
        hippoSubAlv = [1,2,3,4,6,8,9,12,0]
        
    filenamePre = '/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/MultiBlock/' + stain + '/3DVarifold/VarifoldPartFiles/'
    filenamePost = '_FineDown4_jacobW2_' + str(brainNum) + 'only' + savename
    
    selectRegionOnly(filename,filenamePre+'Amygdala'+filenamePost,amy)
    selectRegionOnly(filename,filenamePre+'HATA'+filenamePost,hata)
    selectRegionOnly(filename,filenamePre+'Amygdala+HATA'+filenamePost,amyHata)
    selectRegionOnly(filename,filenamePre+'ERC+Ext'+filenamePost,ercExt)
    selectRegionOnly(filename,filenamePre+'CA1'+filenamePost,ca1)
    selectRegionOnly(filename,filenamePre+'Subiculum'+filenamePost,sub)
    selectRegionOnly(filename,filenamePre+'Hippocampus_All'+filenamePost,hippo)
    selectRegionOnly(filename,filenamePre+'Hippocampus_Sub'+filenamePost,hippoSub)
    selectRegionOnly(filename,filenamePre+'Hippocampus_All+Alv'+filenamePost,hippoAlv)
    selectRegionOnly(filename,filenamePre+'Hippocampus_Sub+Alv'+filenamePost,hippoSubAlv)
    
    
    return

def selectFour(brainNum,stain,mriSuffix,dateNew):
    '''
    Select regions: 
        - amygdala, erc + extension, ca1, subiculum
    Reduce each slice to these regions only so each subset slice can be read and resampled with an isotropic Gaussian
    '''
    savename = '_'+mriSuffix+'_'+dateNew
    if (brainNum == 2):
        five = [2,3,6,7,13]
        sliceNames0 = ['5','6','7','8','9','10','11','12','13','14','15','16','17','18','19']
        sliceNames1 = ['5','6','7','8','9','10','11','12','13','14','15','16','17','18','19']
        sliceNames2 = ['5','6','7','8','9','10','11']
    elif (brainNum == 5):
        five = [14,5,13,1,12]
        sliceNames0 = ['1','5','10','15','20','25','30','35','40','45','50','55','60']
        sliceNames1 = ['1','5','10','15','20','25','30','35','40','45','50','55','60','65']
        sliceNames2 = ['3','7','12','17','22','27','32','37']
        
    filenamePre = '/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/Block1/' + stain + '/3DVarifold/VarifoldPartFiles/'
    filenamePost = '_FineDown4_jacobW2_' + str(brainNum) + 'only' + savename
    
    for sl in sliceNames0:
        fi = filenamePre + 'L' + sl + filenamePost + '.npz'
        fiNew = filenamePre + 'L' + sl + '_FiveRegions' + filenamePost
        selectRegionOnlyInter(fi,fiNew,five)
    
    filenamePre = '/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/Block2/' + stain + '/3DVarifold/VarifoldPartFiles/'
    for sl in sliceNames1:
        fi = filenamePre + 'L' + sl + filenamePost + '.npz'
        fiNew = filenamePre + 'L' + sl + '_FiveRegions' + filenamePost
        selectRegionOnlyInter(fi,fiNew,five)
  
    filenamePre = '/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/Block3/' + stain + '/3DVarifold/VarifoldPartFiles/'
    for sl in sliceNames2:
        fi = filenamePre + 'L' + sl + filenamePost + '.npz'
        fiNew = filenamePre + 'L' + sl + '_FiveRegions' + filenamePost
        selectRegionOnlyInter(fi,fiNew,five)
  
    
    return




def sampleToSurfaces(brainNum,stain,mriSuffix,dateNew,redone):
    '''
    Sample individual particle files to Surfaces using NN
    Surfaces (02/25/22) include:
        amygdala
        ERC + Extension
        CA1
        Subiculum
    '''
    
    savename = '_'+mriSuffix+'_'+dateNew
    fKern = Id
    fKernParams = None
    sKern = nearestNeighbor
    filenamePre = '/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/MultiBlock/' + stain + '/3DVarifold/VarifoldPartFiles/'
    filenamePost = '_FineDown4_jacobW2_' + str(brainNum) + 'only' + savename
    
    sKernParamsA = ['/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/MultiBlock/'+stain+'/3DVarifold/Amygdala_NN_Distances.npy']
    sKernParamsE = ['/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/MultiBlock/'+stain+'/3DVarifold/ERC+Ext_NN_Distances.npy']
    sKernParamsC = ['/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/MultiBlock/'+stain+'/3DVarifold/CA1_NN_Distances.npy']
    sKernParamsS = ['/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/MultiBlock/'+stain+'/3DVarifold/Subiculum_NN_Distances.npy']
    

    
    if (brainNum == 2):
        amyCoords = sp.io.loadmat('/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/3DSegmentations/' + redone + '/Surfaces/amygdala_toMai.mat',appendmat=False,struct_as_record=False)
        amyCoords = np.asarray(amyCoords['YXZ']).astype('float32')
        ercCoords = sp.io.loadmat('/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/3DSegmentations/' + redone + '/Surfaces/ERC+Ext_label01_toMai_0803.mat',appendmat=False,struct_as_record=False)
        ercCoords = np.asarray(ercCoords['YXZ']).astype('float32')
        ca1Coords = sp.io.loadmat('/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/3DSegmentations/' + redone + '/Surfaces/ca1_toMai.mat',appendmat=False,struct_as_record=False)
        ca1Coords = np.asarray(ca1Coords['YXZ']).astype('float32')
        subCoords = sp.io.loadmat('/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/3DSegmentations/' + redone + '/Surfaces/subiculum_toMai.mat',appendmat=False,struct_as_record=False)
        subCoords = np.asarray(subCoords['YXZ']).astype('float32')

    elif (brainNum == 5):
        amyCoords = sp.io.loadmat('/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain5/3DSegmentations/Surfaces/amygdala_toMai_1127.mat',appendmat=False,struct_as_record=False)
        amyCoords = np.asarray(amyCoords['YXZ']).astype('float32')
        ercCoords = sp.io.loadmat('/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain5/3DSegmentations/Surfaces/extension+erc_toMai_1127.mat',appendmat=False,struct_as_record=False)
        ercCoords = np.asarray(ercCoords['YXZ']).astype('float32')
        ca1Coords = sp.io.loadmat('/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain5/3DSegmentations/Surfaces/ca1_toMai_1127.mat',appendmat=False,struct_as_record=False)
        ca1Coords = np.asarray(ca1Coords['YXZ']).astype('float32')
        subCoords = sp.io.loadmat('/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain5/3DSegmentations/Surfaces/subiculum_toMai_1127.mat',appendmat=False,struct_as_record=False)
        subCoords = np.asarray(subCoords['YXZ']).astype('float32')
        
    # Load Varifolds 
    amyfLocT,amyfWeightsT,amyfFeatsT = loadVarifold(filenamePre+'Amygdala'+filenamePost+'.npz')
    print("amygdala feature length")
    print(len(amyfFeatsT))
    for a in amyfFeatsT:
        print(a.shape)
    amyweightsT,amyfeatsT = fineToCoarse(amyCoords,amyfLocT,amyfFeatsT,amyfWeightsT,sKern,sKernParamsA,interToCoarseExVivoHuman11TnonLinear,fKernParams)
    print("amygdala feature length after fine to coarse")
    print(len(amyfeatsT))
    for a in amyfeatsT:
        print(a.shape)
    #np.savez('/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/MultiBlock/' + stain+'/3DVarifold/VarifoldPartFiles/Amygdala_FineDown4_jacobW2_' + str(brainNum) + 'only' + savename + '_Surface_NN.npz',coords=amyCoords,weights=amyweightsT,feats=amyfeatsT)
    feat_dict = dict()
    i = 0
    for f in amyfeatsT:
        feat_dict['feats' + str(i)] = f
        print(f.shape)
        i+=1
    np.savez('/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/MultiBlock/' + stain+'/3DVarifold/VarifoldPartFiles/Amygdala_FineDown4_jacobW2_' + str(brainNum) + 'only' + savename + '_Surface_NN_Dict.npz',coords=amyCoords,weights=amyweightsT,**feat_dict)


    ercfLocT,ercfWeightsT,ercfFeatsT = loadVarifold(filenamePre + 'ERC+Ext'+filenamePost+'.npz')
    print("ERC feature length")
    print(len(ercfFeatsT))
    for e in ercfFeatsT:
        print(e.shape)
    ercweightsT,ercfeatsT = fineToCoarse(ercCoords,ercfLocT,ercfFeatsT,ercfWeightsT,sKern,sKernParamsE,interToCoarseExVivoHuman11TnonLinear,fKernParams)
    print("ERC features after")
    print(len(ercfeatsT))
    for e in ercfeatsT:
        print(e.shape)
    #np.savez('/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/MultiBlock/' + stain+'/3DVarifold/VarifoldPartFiles/ERC+Ext_FineDown4_jacobW2_' + str(brainNum) + 'only' + savename + '_Surface_NN.npz',coords=ercCoords,weights=ercweightsT,feats=ercfeatsT)
    feat_dict = dict()
    i = 0
    for f in ercfeatsT:
        feat_dict['feats' + str(i)] = f
        print(f.shape)
        i+=1
    np.savez('/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/MultiBlock/' + stain+'/3DVarifold/VarifoldPartFiles/ERC+Ext_FineDown4_jacobW2_' + str(brainNum) + 'only' + savename + '_Surface_NN_Dict.npz',coords=ercCoords,weights=ercweightsT,**feat_dict)



    ca1fLocT,ca1fWeightsT,ca1fFeatsT = loadVarifold(filenamePre + 'CA1' + filenamePost+'.npz')
    ca1weightsT,ca1featsT = fineToCoarse(ca1Coords,ca1fLocT,ca1fFeatsT,ca1fWeightsT,sKern,sKernParamsC,interToCoarseExVivoHuman11TnonLinear,fKernParams)
    #np.savez('/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/MultiBlock/' + stain+'/3DVarifold/VarifoldPartFiles/CA1_FineDown4_jacobW2_' + str(brainNum) + 'only' + savename + '_Surface_NN.npz',coords=ca1Coords,weights=ca1weightsT,feats=ca1featsT)
    feat_dict = dict()
    i = 0
    for f in ca1featsT:
        feat_dict['feats' + str(i)] = f
        print(f.shape)
        i+=1
    np.savez('/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/MultiBlock/' + stain+'/3DVarifold/VarifoldPartFiles/CA1_FineDown4_jacobW2_' + str(brainNum) + 'only' + savename + '_Surface_NN_Dict.npz',coords=ca1Coords,weights=ca1weightsT,**feat_dict)


    subfLocT,subfWeightsT,subfFeatsT = loadVarifold(filenamePre + 'Subiculum' + filenamePost+'.npz')
    subweightsT,subfeatsT = fineToCoarse(subCoords,subfLocT,subfFeatsT,subfWeightsT,sKern,sKernParamsS,interToCoarseExVivoHuman11TnonLinear,fKernParams)
    #np.savez('/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/MultiBlock/' + stain+'/3DVarifold/VarifoldPartFiles/Subiculum_FineDown4_jacobW2_' + str(brainNum) + 'only' + savename + '_Surface_NN.npz',coords=subCoords,weights=subweightsT,feats=subfeatsT)
    feat_dict = dict()
    i = 0
    for f in subfeatsT:
        feat_dict['feats' + str(i)] = f
        print(f.shape)
        i+=1
    np.savez('/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/MultiBlock/' + stain+'/3DVarifold/VarifoldPartFiles/Subiculum_FineDown4_jacobW2_' + str(brainNum) + 'only' + savename + '_Surface_NN_Dict.npz',coords=subCoords,weights=subweightsT,**feat_dict)


    return

def writeSurfacesAndLB(brainNum,stain,mriSuffix,dateNew,redone):
    '''
    Sample individual particle files to Surfaces using NN
    Surfaces (02/25/22) include:
        amygdala
        ERC + Extension
        CA1
        Subiculum
    '''
    
    savename = '_'+mriSuffix+'_'+dateNew
    filenamePre = '/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/MultiBlock/' + stain + '/3DVarifold/VarifoldPartFiles/'
    filenamePost = '_FineDown4_jacobW2_' + str(brainNum) + 'only' + savename    
    
    if (brainNum == 2):
        amyF = '/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/3DSegmentations/' + redone + '/Surfaces/amygdala_toMai'
        amyCoords = sp.io.loadmat(amyF+'_LBAll.mat',appendmat=False,struct_as_record=False)
        amyA = np.asarray(amyCoords['A']).astype('float32')
        amyB = np.asarray(amyCoords['B']).astype('float32')
        amyD = np.asarray(amyCoords['D']).astype('float32')
        
        ercF = '/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/3DSegmentations/' + redone + '/Surfaces/ERC+Ext_label01_toMai_0803'
        ercCoords = sp.io.loadmat(ercF+'_LBAll.mat',appendmat=False,struct_as_record=False)
        ercA = np.asarray(ercCoords['A']).astype('float32')
        ercB = np.asarray(ercCoords['B']).astype('float32')
        ercD = np.asarray(ercCoords['D']).astype('float32')
        
        ca1F = '/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/3DSegmentations/' + redone + '/Surfaces/ca1_toMai'
        ca1Coords = sp.io.loadmat(ca1F+'_LBAll.mat',appendmat=False,struct_as_record=False)
        ca1A = np.asarray(ca1Coords['A']).astype('float32')
        ca1B = np.asarray(ca1Coords['B']).astype('float32')
        ca1D = np.asarray(ca1Coords['D']).astype('float32')
        
        subF = '/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/3DSegmentations/' + redone + '/Surfaces/subiculum_toMai'
        subCoords = sp.io.loadmat(subF + '_LBAll.mat',appendmat=False,struct_as_record=False)
        subA = np.asarray(subCoords['A']).astype('float32')
        subB = np.asarray(subCoords['B']).astype('float32')
        subD = np.asarray(subCoords['D']).astype('float32')

    elif (brainNum == 5):
        amyF = '/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain5/3DSegmentations/Surfaces/amygdala_toMai_1127'
        amyCoords = sp.io.loadmat(amyF+'_LBAll.mat',appendmat=False,struct_as_record=False)
        amyA = np.asarray(amyCoords['A']).astype('float32')
        amyB = np.asarray(amyCoords['B']).astype('float32')
        amyD = np.asarray(amyCoords['D']).astype('float32')
        
        ercF = '/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain5/3DSegmentations/Surfaces/extension+erc_toMai_1127'
        ercCoords = sp.io.loadmat(ercF+'_LBAll.mat',appendmat=False,struct_as_record=False)
        ercA = np.asarray(ercCoords['A']).astype('float32')
        ercB = np.asarray(ercCoords['B']).astype('float32')
        ercD = np.asarray(ercCoords['D']).astype('float32')

        ca1F = '/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain5/3DSegmentations/Surfaces/ca1_toMai_1127'
        ca1Coords = sp.io.loadmat(ca1F + '_LBAll.mat',appendmat=False,struct_as_record=False)
        ca1A = np.asarray(ca1Coords['A']).astype('float32')
        ca1B = np.asarray(ca1Coords['B']).astype('float32')
        ca1D = np.asarray(ca1Coords['D']).astype('float32')
        
        subF = '/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain5/3DSegmentations/Surfaces/subiculum_toMai_1127'
        subCoords = sp.io.loadmat(subF+'_LBAll.mat',appendmat=False,struct_as_record=False)
        subA = np.asarray(subCoords['A']).astype('float32')
        subB = np.asarray(subCoords['B']).astype('float32')
        subD = np.asarray(subCoords['D']).astype('float32')

    subTot = '/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/MultiBlock/' + stain+'/3DVarifold/VarifoldPartFiles/Subiculum_FineDown4_jacobW2_' + str(brainNum) + 'only' + savename + '_Surface_NN_Dict'
    amyTot = '/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/MultiBlock/' + stain+'/3DVarifold/VarifoldPartFiles/Amygdala_FineDown4_jacobW2_' + str(brainNum) + 'only' + savename + '_Surface_NN_Dict'
    ca1Tot = '/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/MultiBlock/' + stain+'/3DVarifold/VarifoldPartFiles/CA1_FineDown4_jacobW2_' + str(brainNum) + 'only' + savename + '_Surface_NN_Dict'
    ercTot = '/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/MultiBlock/' + stain+'/3DVarifold/VarifoldPartFiles/ERC+Ext_FineDown4_jacobW2_' + str(brainNum) + 'only' + savename + '_Surface_NN_Dict'
    vf.writeRegionVarifoldToVTK(amyTot+'.npz',amyF,-1,"ERC",amyTot+'_LBAll_k2.vtk',amyB,amyD,amyA,2)
    vf.writeRegionVarifoldToVTK(ercTot+'.npz',ercF,-1,"ERC",ercTot+'_LBAll_k2.vtk',ercB,ercD,ercA,2)
    vf.writeRegionVarifoldToVTK(ca1Tot+'.npz',ca1F,-1,"ERC",ca1Tot+'_LBAll_k2.vtk',ca1B,ca1D,ca1A,2)
    vf.writeRegionVarifoldToVTK(subTot+'.npz',subF,-1,"ERC",subTot+'_LBAll_k2.vtk',subB,subD,subA,2)
    
    return
    
def resampleIsotropicGauss(brainNum,stain,mriSuffix,dateNew):
    '''
    brainNum is 1 based here
    '''
    sigma = [0.1,0.1,0.1] # default is isotropic Gaussian
    intdenom = 0.4 # resample at intervals of 0.25 mm
    numDown = 4
    numLab = 5 # amygdala, ERC, extension, CA1, subiculum (no background here)
    fromFine = True
    gauss=True
    dims = None
    byLabel=True
    toVivo = True
    minX = [-30,14,-3]
    maxX = [-17,31,22]
    minX = [-26,7,0]
    maxX = [5,35,42]
    suff = '_' + mriSuffix + '_' + dateNew
    weightBy = '_jacobW2_' + str(brainNum) + 'only' # assume 2D weights
    subsetSuff = weightBy + suff
    
    if (stain == 'Tau'):
        stainNum = 2
    elif (stain == 'Amyloid'):
        stainNum = 0
    
    # get Coarse grid based on individual 2D to 3D points
    saveName = '/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/MultiBlock/' + stain + '/3DVarifold/CoarseGrid_1-4mm' + suff + '.npz'
    
    regionString = '_FiveRegions'
    
    getCoarseGrid(brainNum,[0,1,2],stainNum,sigma,intdenom,dims,minX,maxX,saveName,subsetSuff)
    scaleGauss(saveName,brainNum,[0,1,2],sigma,numLab,stainNum,numDown,intdenom,gauss,byLabel,subsetSuff,regionString,toVivo) # coarsest grid every 4 mm
    
    return


    
    