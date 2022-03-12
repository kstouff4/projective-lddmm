################################################################
# Author: Katie Stouffer, 8/22/21

# Tools for constructing, manipulating, sampling, etc. varifolds. These tools come mostly from the General_Varifold and Varifold framework python notebooks.

################################################################

import sys
sys.path.append('../')
sys.path.append('../SurfaceTools/')

import deformSegmentations3D as ds
import vtkFunctions as vt

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

from matplotlib.lines import Line2D

###############################################################
# Functions for Loading Varifolds

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


############################################################
# Functions for Altering Varifolds
def addFeatures(varFile,feats,savename):
    fpLoc,fpW,fpF = loadVarifold(varFile)
    fpF.append(feats)
    feats_dict = dict()
    for i in range(len(fpF)):
        feats_dict['feats'+str(i)] = fpF[i]
    np.savez(savename,coords=fpLoc,weights=fpW,**feats_dict)
    return

def addLabelsFromVolToMaiSpace(varFile,volLabFile,savename,rigidMat):
    '''
    Assume volLabFile generated from deformSegmentations3D so contains allVolLabs, maskedVolLabs, and coords (in MRI space)
    Assume d0=d1=d2 = 0.125 for coordinates
    '''
    
    fpLoc,fpW,fpF = loadVarifold(varFile)
    
    # move rigid backwards from mai to MRI
    paramsFirst = sp.io.loadmat(rigidMat,appendmat=False,struct_as_record=False)
    A = np.asarray(paramsFirst['A'])
    # transpose 
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
       
    Ainv = np.linalg.inv(Aret)
    coordsRet = np.zeros_like(fpLoc)
    coordsRet[...,0] = Ainv[0,0]*fpLoc[...,0] + Ainv[0,1]*fpLoc[...,1] + Ainv[0,2]*fpLoc[...,2] + Ainv[0,3]
    coordsRet[...,1] = Ainv[1,0]*fpLoc[...,0] + Ainv[1,1]*fpLoc[...,1] + Ainv[1,2]*fpLoc[...,2] + Ainv[1,3]
    coordsRet[...,2] = Ainv[2,0]*fpLoc[...,0] + Ainv[2,1]*fpLoc[...,1] + Ainv[2,2]*fpLoc[...,2] + Ainv[2,3]
    
    vol = np.load(volLabFile)
    coords = vol['coords'].astype('float32')
    maskLabs = vol['maskedVolLabs'].astype('float32')
    fpNewLab = ds.applyFunction(coordsRet,maskLabs,np.asarray([coords[0,0,0,0],coords[1,0,0,0]]),np.asarray([coords[0,0,0,1],coords[0,1,0,1]]),np.asarray([coords[0,0,0,2],coords[0,0,1,2]]),NN=2)
    
    addFeatures(varFile,fpNewLab,savename)
    return


###########################################################
# Functions for Analyzing Varifolds

def makeBarFromVarifold(brainNum,stainString,exp,weightBy,lab2D=False,normalize=False,factor=1.0,ylim=None):
    # Make color dictionary consistent
    labelsD = ['other','HATA','alveus','amygdala','ca1','ca2','ca3','ERC-extension','ERC','DG-granular','DG-hilus','DG-molecular','parasubiculum','presubiculum','subiculum','endfolial']
    colorsD = ['black','olive','red','green','lime','cyan','yellow','thistle','plum','orange','cornsilk','goldenrod','slateblue', 'purple', 'fuchsia','navy']
    cDict = {labelsD[i].lower():colorsD[i] for i in range(len(labelsD))}
    
    
    if (brainNum == 2 and lab2D):
        saveNamePost = '_FineDown' + str(exp) + '_2DLabels_jacobW' + str(weightBy) + '.npz'
        labels = ['alveus','amygdala','CA1','CA2','CA3','erc','dentate_granular','hilus','dentate_molecular','parasubiculum','presubiculum','subiculum','HATA','TEC','other']
        labels = ['HATA','Alveus','Amygdala','CA1','CA2','CA3','ERC-extension','ERC','DG-granular','DG-hilus','DG-molecular','Parasubiculum','Presubiculum','Subiculum','other']
        orderingR = [14,13,12,1,0,5,9,10,11,2,3,4,8,6,7]
        ordering = [7,6,8,4,3,2,11,10,9,5,0,1,12,13,14]
        sliceNames = [5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]
        amy = 2
    elif (brainNum == 2):
        saveNamePost = '_FineDown' + str(exp) + '_jacobW' + str(weightBy) + '.npz'
        labels = ['alveus','CA1','CA2','CA3','erc','DG_granular','DG_hilus','DG_molecular','parasubiculum','presubiculum','subiculum','tec','amygdala','other']
        labels = ['HATA','Alveus','Amygdala','CA1','CA2','CA3','ERC-extension','ERC','DG-granular','DG-hilus','DG-molecular','Parasubiculum','Presubiculum','Subiculum','other']
        ordering = [12,11,4,8,9,10,1,2,3,7,5,6,0,13]
        orderingR = [13,0,6,5,7,3,2,1,10,9,8,4,11,12]
        ordering = [2,0,6,7,11,12,13,3,4,5,10,8,9,1,14]
        orderingR = [14,1,9,8,10,5,4,3,13,12,11,7,6,0,2]
        sliceNames = [5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]
        amy = 2
        erc = 7
    elif (brainNum == 5):
        labels = ['molecular','granular','ca2','ca1','erc','ca3','srlm+VHS','subiculum','alveus','fasciolar','hilus','endfolial','presubiculum','parasubiculum','HATA','other']
        labels = ['alveus','ca1','ca2','ca3','endfolial','erc','granular','HATA','hilus','molecular','parasubiculum','presubiculum','subiculum','other']
        labels = ['alveus','ca1','ca2','ca3','endfolial','erc','DG-granular','HATA','DG-hilus','DG-molecular','parasubiculum','presubiculum','subiculum','ERC-extension','amygdala','other']
        amy = 14
        saveNamePost = '_FineDown' + str(exp) + '_jacobW' + str(weightBy) + '.npz'
        ordering = [14,7,13,5,10,11,12,1,2,3,4,9,6,8,0,15]
        orderingR = [15,0,8,6,9,4,3,2,1,12,11,10,5,13,7,14]
        sliceNames = [1,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75]
    elif (brainNum == 3):
        labels = ['alveus','ca2','ca3','endfolial','granular','hilus','molecular','ca1','subiculum','erc','HATA','presubiculum','parasubiculum','amygdala','other']
        saveNamePost = '_FineDown' + str(exp) + '_jacobW' + str(weightBy) + '.npz'
        amy = 13
        orderingR = [14,13,10,0,9,12,11,8,7,1,2,3,6,4,5]
        ordering = [5,4,6,3,2,1,7,8,11,12,9,0,10,13,14]
        sliceNames = [1,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75]

    print("save name is ")
    print(saveNamePost)
    # calculate total tau and total area 
    numLab = len(labels)
    tauPerLab = np.zeros((numLab,1))
    tauPerLabNum = np.zeros((numLab,1))
    tauPerLabDen = np.zeros((numLab,1))
    for blockNum in range(3):
        high = 15
        if (brainNum == 5 and blockNum == 2):
            high = 8
            sliceNames = [3,7,12,17,22,27,32,37]
        elif (brainNum == 3 and ((blockNum == 0) or (blockNum == 2))):
            high = 13
        elif (brainNum == 3 and blockNum == 1):
            high = 12
        elif (brainNum == 3 and blockNum == 3):
            high = 11
        for sl in range(high):
            locationNum = sliceNames[sl]
            savename = '/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/Block' + str(blockNum+1) + '/' + stainString + '/3DVarifold/VarifoldPartFiles/L' + str(locationNum) + saveNamePost
            print(savename)
            if (not path.exists(savename)):
                continue
            loc,weights,feats = loadVarifold(savename)
            print("path exists")
            for l in range(numLab):
                tauPerLabNum[l] += np.sum(weights[...,None]*feats[1][...,l])
                tauPerLabDen[l] += np.sum(weights[...,None]*feats[2][...,l])
    
    tauPerLabDenR = 1.0/tauPerLabDen
    tauPerLabDenR[tauPerLabDenR == np.infty] = 0
    tauPerLab = tauPerLabNum*tauPerLabDenR
    f,ax = plt.subplots(figsize=(6,10),constrained_layout=True)
    for lab in range(numLab):
        labI = int(ordering[lab])
        ax.barh(numLab-lab,tauPerLab[labI]*factor,color=cDict[labels[labI].lower()])
    
    ax.set_ylabel('Region')
    ax.set_yticks(list(np.arange(len(orderingR))+1))
    ax.set_yticklabels([labels[t] for t in orderingR])
    ax.set_xlabel('Tau Tangles / mm$^2$')
    ax.set_title("Average Density of " + stainString + " Particles By Region Across Blocks")
    if (ylim is not None):
        ax.set_xlim(0,ylim)
    #ax.set_xlim(0,150)
    f.canvas.draw()
    f.savefig('/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/MultiBlock/' + stainString + '/BarAndLine/barGraphFromVarifold' + str(weightBy) + '.png',dpi=300)
    if (normalize):
        f,ax = plt.subplots(figsize=(6,10),constrained_layout=True)
        for lab in range(numLab):
            amy = np.argmax(tauPerLab)
            labI = int(ordering[lab])
            ax.barh(numLab-lab,tauPerLab[labI]/tauPerLab[amy],color=cDict[labels[labI].lower()])

            ax.set_ylabel('Region')
            ax.set_yticks(list(np.arange(len(orderingR))+1))
            ax.set_yticklabels([labels[t] for t in orderingR])
            ax.set_xlabel('Number of particles / mm$^2$, normalized to ' + str(labels[amy]) + ' density')
            ax.set_title("Average Density of " + stainString + " Particles By Region Across Blocks")
            ax.set_xlim(0,1)
            f.canvas.draw()
            f.savefig('/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/MultiBlock/' + stainString + '/BarAndLine/barGraphFromVarifold' + str(weightBy) + '_normalizedToERC.png',dpi=300)

    return

def makeDoubleBarFromVarifold(brainNums,stainString,exp,weightBys,lab2D=False,normalize=False,factors=[1.0,1.0],ylim=None):
    # Make color dictionary consistent
    labelsD = ['other','HATA','alveus','amygdala','ca1','ca2','ca3','ERC-extension','ERC','DG-granular','DG-hilus','DG-molecular','parasubiculum','presubiculum','subiculum','endfolial']
    colorsD = ['black','olive','red','green','lime','cyan','yellow','thistle','plum','orange','cornsilk','goldenrod','slateblue', 'purple', 'fuchsia','navy']
    cDict = {labelsD[i].lower():colorsD[i] for i in range(len(labelsD))}
    
    b = 0
    offs = [-0.2,-0.8]
    edges = ['black','gray']
    leg = ['Brain Sample 1','Brain Sample 2']
    custom_lines = [Line2D([0], [0], color='black', lw=4),
                Line2D([0], [0], color='gray', lw=4)]
    f,ax = plt.subplots(figsize=(10,6),constrained_layout=True)
    for brainNum in brainNums:
        factor = factors[b]
        weightBy = weightBys[b]
        if (brainNum == 2 and lab2D):
            saveNamePost = '_FineDown' + str(exp) + '_2DLabels_jacobW' + str(weightBy) + '.npz'
            labels = ['alveus','amygdala','CA1','CA2','CA3','erc','dentate_granular','hilus','dentate_molecular','parasubiculum','presubiculum','subiculum','HATA','TEC','other']
            labels = ['HATA','Alveus','Amygdala','CA1','CA2','CA3','ERC-extension','ERC','DG-granular','DG-hilus','DG-molecular','Parasubiculum','Presubiculum','Subiculum','other']
            orderingR = [14,13,12,1,0,5,9,10,11,2,3,4,8,6,7]
            ordering = [7,6,8,4,3,2,11,10,9,5,0,1,12,13,14]
            sliceNames = [5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]
            amy = 2
        elif (brainNum == 2):
            saveNamePost = '_FineDown' + str(exp) + '_jacobW' + str(weightBy) + '.npz'
            labels = ['alveus','CA1','CA2','CA3','erc','DG_granular','DG_hilus','DG_molecular','parasubiculum','presubiculum','subiculum','tec','amygdala','other']
            labels = ['HATA','Alveus','Amygdala','CA1','CA2','CA3','ERC-extension','ERC','DG-granular','DG-hilus','DG-molecular','Parasubiculum','Presubiculum','Subiculum','other']
            ordering = [12,11,4,8,9,10,1,2,3,7,5,6,0,13]
            orderingR = [13,0,6,5,7,3,2,1,10,9,8,4,11,12]
            ordering = [2,0,6,7,11,12,13,3,4,5,10,8,9,1,14]
            orderingR = [14,1,9,8,10,5,4,3,13,12,11,7,6,0,2]
            sliceNames = [5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]
            amy = 2
            erc = 7
        elif (brainNum == 5):
            labels = ['molecular','granular','ca2','ca1','erc','ca3','srlm+VHS','subiculum','alveus','fasciolar','hilus','endfolial','presubiculum','parasubiculum','HATA','other']
            labels = ['alveus','ca1','ca2','ca3','endfolial','erc','granular','HATA','hilus','molecular','parasubiculum','presubiculum','subiculum','other']
            labels = ['alveus','ca1','ca2','ca3','endfolial','erc','DG-granular','HATA','DG-hilus','DG-molecular','parasubiculum','presubiculum','subiculum','ERC-extension','amygdala','other']
            amy = 14
            saveNamePost = '_FineDown' + str(exp) + '_jacobW' + str(weightBy) + '.npz'
            ordering = [14,7,13,5,10,11,12,1,2,3,4,9,6,8,0,15]
            orderingR = [15,0,8,6,9,4,3,2,1,12,11,10,5,13,7,14]
            sliceNames = [1,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75]
        elif (brainNum == 3):
            labels = ['alveus','ca2','ca3','endfolial','granular','hilus','molecular','ca1','subiculum','erc','HATA','presubiculum','parasubiculum','amygdala','other']
            saveNamePost = '_FineDown' + str(exp) + '_jacobW' + str(weightBy) + '.npz'
            amy = 13
            orderingR = [14,13,10,0,9,12,11,8,7,1,2,3,6,4,5]
            ordering = [5,4,6,3,2,1,7,8,11,12,9,0,10,13,14]
            sliceNames = [1,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75]

        print("save name is ")
        print(saveNamePost)
        # calculate total tau and total area 
        numLab = len(labels)
        tauPerLab = np.zeros((numLab,1))
        tauPerLabNum = np.zeros((numLab,1))
        tauPerLabDen = np.zeros((numLab,1))
        for blockNum in range(3):
            high = 15
            if (brainNum == 5 and blockNum == 2):
                high = 8
                sliceNames = [3,7,12,17,22,27,32,37]
            elif (brainNum == 3 and ((blockNum == 0) or (blockNum == 2))):
                high = 13
            elif (brainNum == 3 and blockNum == 1):
                high = 12
            elif (brainNum == 3 and blockNum == 3):
                high = 11
            for sl in range(high):
                locationNum = sliceNames[sl]
                savename = '/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNum) + '/Block' + str(blockNum+1) + '/' + stainString + '/3DVarifold/VarifoldPartFiles/L' + str(locationNum) + saveNamePost
                print(savename)
                if (not path.exists(savename)):
                    continue
                loc,weights,feats = loadVarifold(savename)
                print("path exists")
                for l in range(numLab):
                    tauPerLabNum[l] += np.sum(weights[...,None]*feats[1][...,l])
                    tauPerLabDen[l] += np.sum(weights[...,None]*feats[2][...,l])

        tauPerLabDenR = 1.0/tauPerLabDen
        tauPerLabDenR[tauPerLabDenR == np.infty] = 0
        tauPerLab = tauPerLabNum*tauPerLabDenR
        for lab in range(numLab):
            labI = int(orderingR[lab])
            if (labels[labI].lower() == "endfolial"):
                offs[b] = offs[b]+1
                print("changing offset")
                continue
            ax.bar(numLab-lab + offs[b],tauPerLab[labI]*factor,width=0.4,color=cDict[labels[labI].lower()],edgecolor=edges[b],hatch="/")
        b += 1

    labels = ['HATA','Alveus','Amygdala','CA1','CA2','CA3','ERC-extension','ERC','DG-granular','DG-hilus','DG-molecular','Parasubiculum','Presubiculum','Subiculum','other']
    orderingR = [14,1,9,8,10,5,4,3,13,12,11,7,6,0,2]
    ordering = [2,0,6,7,11,12,13,3,4,5,10,8,9,1,14]
    ax.set_xlabel('Region')
    ax.set_xticks(list(np.arange(len(ordering))+1))
    ax.set_xticklabels([labels[t] for t in ordering])
    for tick in ax.get_xticklabels():
        tick.set_rotation(45)
    ax.set_ylabel('Tau Tangles / mm$^2$')
    ax.set_title("Average Density of " + stainString + " Particles By Region Across Blocks")
    if (ylim is not None):
        ax.set_ylim(0,ylim)
        #ax.set_xlim(0,150)
    ax.legend(custom_lines, ['Brain Sample 1', 'Brain Sample 2'])
    f.canvas.draw()
    f.savefig('/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain' + str(brainNums[0]) + '/MultiBlock/' + stainString + '/BarAndLine/barGraphFromVarifold' + str(weightBy) + 'horiz.png',dpi=300)


    return




def makeBarFromVarifoldLabel(varFile,fVal,fCat,catNames,catColors,savename,fVal2 = -1, fCat2 = -1,fValD = -1, fValD2 = -1,weight=None):
    '''
    count up fVal in all of the varifold particles per fCat
    if weight = True, then count is fVal*weight / fCat*weight total
    
    Example: fVal = 1,5 = tau/ERC in brain 2
    
    To get standard bargraph run General Varifold
    '''
    fpLoc,fpW,fpF = loadVarifold(varFile)
    print("shapes of loc,w,f")
    print(fpLoc.shape)
    print(fpW.shape)
    print(len(fpF))
    
    fCount = fpF[fVal]
    fCateg = fpF[fCat]
    if (fValD >= 0):
        fDenom = fpF[fValD]
        if (fValD2 >= 0):
            fDenom = fDenom[...,fValD2]
    if (fVal2 >= 0):
        fCount = fCount[...,fVal2] # example = tau / erc 
    if (fCat2 >= 0):
        fCateg = fCateg[...,fCat2] # example = area ERC (will be fraction of area ERC so weight will be 
    countsPerCatN = np.zeros((len(catNames),1))
    countsPerCatD = np.zeros((len(catNames),1))
    fCount = np.ravel(fCount)
    fDenom = np.ravel(fDenom) # assume arrays are of dimension 1 for second dimension
    fCateg = np.ravel(fCateg)
    print("unique categories")
    print(np.unique(fCateg))
    if (weight is not None):
        if (weight[0] == -1):
            fpW = np.ravel(fpW)
            fCount = fCount*fpW
            fDenom = fDenom*fpW
        else:
            weightToMul = fpF[weight[0]]
            if len(weight) > 1:
                weightToMul = weightToMul[...,weight[1]]
            wRav = np.ravel(weightToMul)
            fCount = fCount*wRav*np.ravel(fpW)
            fDenom = fDenom*wRav*np.ravel(fpW)
    for c in range(len(catNames)):
        countsPerCatN[c] = np.sum((fCateg == c)*fCount) # i.e. tau tangles in ERC (in area c)
        countsPerCatD[c] = np.sum((fCateg == c)*fDenom) # i.e. area of ERC (in area c)
    print('Before Divide')
    print(countsPerCatN)
    print(countsPerCatD)
    rat = countsPerCatN/countsPerCatD # assume no 0 in denom
    f,ax = plt.subplots(figsize=(6,10),constrained_layout=True)
    totalC = len(catNames)
    for c in range(totalC):
        ax.barh(c,rat[c],color=catColors[c])
    
    ax.set_ylabel('Category')
    ax.set_yticks(list(np.arange(totalC)))
    ax.set_yticklabels(catNames)
    ax.set_xlabel('Number of particles / mm$^2$')
    ax.set_title("Average Density of Particles By Region Across Blocks")
    ax.set_xlim(0,215)
    f.canvas.draw()
    f.savefig(savename,dpi=300)
    return

def plotSlidingWindow(varFiles,width,savename,axis=2,inc=1,minCoord=100000,maxCoord=-100000,featsN = 0):
    '''
    Plot total tau along a particular axis by summing up between window centered at each point (increment 1)
    where you have width = width on eitehr side of center
    
    '''
    if (minCoord > 10000):
        minCoord = 100000
        maxCoord = -100000
        for vf in varFiles:
            loc,weights,feats = loadVarifold(vf)
            if (np.min(loc[...,axis]) < minCoord):
                minCoord = np.min(loc[...,axis])
            if (np.max(loc[...,axis]) > maxCoord):
                maxCoord = np.max(loc[...,axis])
    minCoord = np.floor(minCoord)
    maxCoord = np.ceil(maxCoord)
    
    inds = np.arange(minCoord + width,maxCoord-width+1,inc)
    indVals = np.zeros_like(inds)
    indArea = np.zeros_like(inds)
    for vf in varFiles:
        loc,weights,feats = loadVarifold(vf)
        print(weights.shape)
        print(loc.shape)
        print(feats[featsN].shape)
        for i in range(len(inds)):
            indVals[i] += np.sum(np.squeeze(feats[featsN])*np.squeeze(weights)*np.squeeze(loc[...,axis] >= inds[i]-width)*np.squeeze(loc[...,axis] < inds[i] + width)) # total tau
            indArea[i] += np.sum(np.squeeze(weights)*np.squeeze(loc[...,axis] >= inds[i]-width)*np.squeeze(loc[...,axis] < inds[i] + width))
    
    recipA = 1.0/indArea
    recipA[recipA == np.infty] = 0
    rat = recipA*indVals
    # plot 
    f,ax = plt.subplots()
    ax.plot(rat,inds) # flipped so is vertical
    ax.set_xlabel('Average Density (per mm^2) of Tau Tangles in ' + str(2*width) + ' mm along Z axis')
    ax.set_ylabel('Center of ' + str(2*width) + ' mm window')
    f.savefig(savename)
    return inds,rat,indArea

def plotTrends(savename,xcoords,vals,areas,legLab=["Amygdala","ERC","Subiculum","CA1"],title="Tau Tangle Density Rostral to Caudal",xlab="Center of 2 mm window (mm, Mai Z Axis)",ylab="Tau Tangles / mm^2",cols=['green','plum','fuchsia','lime']):
    f,ax = plt.subplots(figsize=(10,6))
    ax.set_title(title)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    X = np.ones((xcoords.shape[0],2))
    X[...,1] = xcoords # n x 2
    XT = X.T # 2 x n
    widthC = [0.8/2, 0.8/4, -0.8/4, -0.8/2]
    for i in range(len(vals)):
        v = vals[i]
        print(v.shape)
        lab = legLab[i]
        W = np.zeros((xcoords.shape[0],xcoords.shape[0]))
        a = areas[i]
        xless = xcoords[a > 0]
        numX = xless.shape[0]
        X = np.ones((numX,2))
        X[...,1] = xless # n x 2
        XT = X.T # 2 x n
        vless = v[a > 0]

        #for j in range(xcoords.shape[0]):
        #    W[j,j] = int(a[j] > 0)
        slope = np.linalg.inv(XT@X)@XT@vless #2 x 2 * 2 x 1 = 2 x 1
        if (len(vals) == 1):
            ax.bar(xcoords,v,width=1,color=cols[i],label=lab + ", slope is " + str(np.round(slope[1],4)),alpha=0.5)
        else:
            ax.bar(xcoords-widthC[i],v,width=0.4,color=cols[i],label=lab + ", slope is " + str(np.round(slope[1],4)),alpha=0.5)
        ax.plot(xless,slope[1]*xless + slope[0],linewidth=4,linestyle='--',color=cols[i])
        #ax.text(0.95,0.9,f"Slope is {slope[1]:.4f}",horizontalalignment='center',verticalalignment='center',transform=ax.transAxes,fontsize=8)

    ax.legend()
    f.savefig(savename,dpi=300)
    return

def plotTrendsSeaborn(savename,xcoords,vals,areas,legLab=["Amygdala","ERC","Subiculum","CA1"],title="Tau Tangle Density Rostral to Caudal",xlab="Center of 2 mm window (mm, Mai Z Axis)",ylab="Tau Tangles / mm$^2$",cols=['green','plum','fuchsia','lime'],trend=True,factor=1,xlim=None):
    boot = 1000
    '''      
    datanp = np.zeros((xcoords.shape[0],len(vals)+1))
    datanp[:,0] = xcoords
    for i in range(len(vals)):
        datanp[:,i+1] = vals[i]
    df = pd.DataFrame(data=datanp,columns = legLab.insert(0,'Center of 2 mm window (mm, Mai Z Axis)'))
    '''
    f,ax = plt.subplots(figsize=(10,6))
    ax.set_title(title)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    X = np.ones((xcoords.shape[0],2))
    X[...,1] = xcoords # n x 2
    XT = X.T # 2 x n
    for i in range(len(vals)):
        v = vals[i]*factor
        print(v.shape)
        lab = legLab[i]
        W = np.zeros((xcoords.shape[0],xcoords.shape[0]))
        a = areas[i]
        xless = xcoords[a > 0]
        numX = xless.shape[0]
        X = np.ones((numX,2))
        X[...,1] = xless # n x 2
        XT = X.T # 2 x n
        vless = v[a > 0]

        #for j in range(xcoords.shape[0]):
        #    W[j,j] = int(a[j] > 0)
        if (trend):
            slope = np.linalg.inv(XT@X)@XT@vless #2 x 2 * 2 x 1 = 2 x 1
            ax.bar(xcoords,v,width=1,color=cols[i],label=lab + " Trend Line (slope " + str(np.round(slope[1],4)) + ")",alpha=0.4)
            ax.plot(xless,slope[1]*xless + slope[0],linewidth=4,color=cols[i])
            # bootstrap
            c = xless.shape[0]
            xcoordsB = np.zeros((c,boot))
            valsB = np.zeros((vless.shape[0],boot))
            slopes = np.zeros((boot,1))
            for b in range(boot):
                r = random.choices(np.arange(c),k=c)
                xcoordsB[:,b] = xless[r]
                X = np.ones((numX,2))
                X[...,1] = xcoordsB[:,b]# n x 2
                XT = X.T 
                slope = np.linalg.inv(XT@X)@XT@vless[r] #2 x 2 * 2 x 1 = 2 x 1
                slopes[b] = slope[1]
                xlessO = np.sort(xless[r])
                xcoordsB[:,b] = xlessO
                valsB[:,b] = slope[1]*xless + slope[0]
            lowS = np.quantile(slopes,0.025)
            highS = np.quantile(slopes,0.975)
            lowV = np.quantile(valsB,0.025,axis=-1)
            highV = np.quantile(valsB,0.975,axis=-1)
            ax.plot(xless,lowV,linewidth=2,linestyle=':',color='red',label='95% Confidence Interval\nSlopes: [' + str(np.round(lowS,4)) + ', ' + str(np.round(highS,4)) + ']')
            ax.plot(xless,highV,linewidth=2,linestyle=':',color='red')
        else:
            ax.bar(xcoords,v,width=1,color=cols[i],label=lab,alpha=0.4)
        if (xlim is not None):
            ax.set_ylim([0,xlim])
        #ax.text(0.95,0.9,f"Slope is {slope[1]:.4f}",horizontalalignment='center',verticalalignment='center',transform=ax.transAxes,fontsize=8)

    ax.legend()
    f.savefig(savename,dpi=300)
    return

def plotMaiCoronalTau(varFiles,width,savename,axis=2,center=10,featsN=0,maxVal=150):
    '''
    Plot total tau along a particular axis by summing up between window centered at each point (increment 1)
    where you have width = width on eitehr side of center
    
    '''
     
    indVals = []
    indX = []
    indY = []
    indArea = []
    indValsNum = 0
    indXNum = 0
    indYNum = 0
    indAreaNum = 0
    for vf in varFiles:
        loc,weights,feats = loadVarifold(vf)
        featsToAdd = np.squeeze(weights[(loc[...,axis] >= center-width)*(loc[...,axis] < center + width)])*np.squeeze(feats[featsN][(loc[...,axis] >= center-width)*(loc[...,axis] < center + width)]) # total tau
        xToAdd = np.squeeze(loc[(loc[...,axis] >= center-width)*(loc[...,axis] < center + width),1])
        yToAdd = np.squeeze(loc[(loc[...,axis] >= center-width)*(loc[...,axis] < center + width),0])
        areaToAdd = np.squeeze(weights[(loc[...,axis] >= center-width)*(loc[...,axis] < center + width)])
        indVals.append(featsToAdd)
        indX.append(xToAdd)
        indY.append(yToAdd)
        indArea.append(areaToAdd)
        indValsNum += featsToAdd.shape[0]
        indXNum += xToAdd.shape[0]
        indYNum += yToAdd.shape[0]
        indAreaNum += areaToAdd.shape[0]
        
    indValsA = np.zeros((indValsNum,1))
    indXA = np.zeros((indXNum,1))
    indYA = np.zeros((indYNum,1))
    indAreaA = np.zeros((indAreaNum,1))
    
    i = 0
    for f in range(len(indVals)):
        iv = indVals[f]
        s = iv.shape[0]
        indValsA[i:i+s] = iv[...,None]
        
        ix = indX[f]
        indXA[i:i+s] = ix[...,None]
        
        iy = indY[f]
        indYA[i:i+s] = iy[...,None]
        
        ia = indArea[f]
        indAreaA[i:i+s] = ia[...,None]
        i = i+s
   

    print(indValsA.shape)
    print(indXA.shape)
    recipA = 1.0/indAreaA
    recipA[recipA == np.infty] = 0
    rat = recipA*indValsA
    print(rat.shape)
        
    f,ax = plt.subplots()
    ax.set_title("Mai Slice Centered at " + str(center) + " mm of width " + str(width))
    im = ax.scatter(indXA,indYA,s=indAreaA*10,c=rat,cmap='jet',vmax=np.quantile(rat,0.96),vmin=0)
    c = f.colorbar(im,ax=ax[0])
    c.set_label('Tau Tangles / mm$^2$')
    #ax[1].set_title("Tau Only Mai Slice " + str(center))
    #im = ax[1].scatter(indXA,indYA,s=indAreaA*10,c=rat,cmap='jet',vmax=np.quantile(rat,0.96),vmin=1)
    #c = f.colorbar(im,ax=ax[1])
    #c.set_label('Tau Tangles/mm$^2$')
    f.savefig(savename)
    
    return indValsA,indAreaA,indXA,indYA

        

################################################################
# Functions for Converting Varifolds to Other Formats

def writeRegionVarifoldToVTK(varFile,pointMatFile,regNum,regName,savename,B=None,D=None,A=None,k=1):
    
    # read in YXZ and polys data (assume polys is 1 based and needs to be revised)
    params = sp.io.loadmat(pointMatFile)
    polys = np.asarray(params['polys'])
    YXZ = np.asarray(params['YXZ'])
    if (np.min(polys) == 1):
        polys[:,1:] = polys[:,1:]-1
    
    loc,weights,feats = loadVarifold(varFile)
    
    featsVTK = []
    featsVTKNames = []
    
    # smoothe with LB basis
    if (B is None):
        LB = np.eye(feats[0].shape[0])
        A = np.eye(feats[0].shape[0])
    else:
        #LB = B@B.T # puts it into a new basis
        if D is None:
            Dmat = np.eye(B.shape[1])
        else:
            Dmat = (D)*np.eye(B.shape[1]) # e^-\lambda^2 along the diagonals
            verts = B.shape[1]
            print("B and D check")
            print(B[0:3,0:3])
            print(D[0:3,0:3])
            DmatNew = np.zeros_like(Dmat)
            print("shape of Dmat")
            print(Dmat.shape)
            print("A shape")
            print(A.shape)
            for r in range(verts):
                for c in range(verts):
                    DmatNew[r,c] = 1.0/(1.0-k*Dmat[c,c]*A[r,r])
            #part = np.sqrt(1.0/np.sum(Dmat*Dmat)) # sum of all components should be 1
            #Dmat = Dmat*part
        #Bnorm = (1.0/np.sum(B*B,axis=0))*B
        print("Coefficients")
        print(np.min(DmatNew))
        print(np.max(DmatNew))
        check = np.unique(B.T@A@B)
        print("Check for orthogonality according to A norm: ")
        print(check)
        print(np.min(B.T@A@B))
        print(np.max(B.T@A@B))
        Bscale = DmatNew*B
        LB = Bscale@B.T@A
        LBin = np.linalg.inv(LB)
        LBId = B@(np.eye(B.shape[1]))@np.linalg.inv(B)
    # add Weights of Area 
    if (regNum >= 0):
        featsVTK.append(LB@feats[2][:,regNum])
        featsVTKNames.append('FRACTION_OF_AREA_' + str(regName))
    
        # add Tau in Area 
        featsVTK.append(LB@feats[3][:,regNum])
        featsVTKNames.append('TAU_TANGLES_PER_MMSQ_' + str(regName))
    
    # add Total Weight 
    #b = vt.computeSurfaceIntegral(weights,YXZ,polys)
    #a = vt.computeSurfaceIntegral(LB@weights,YXZ,polys)
    #b = np.sum(A@weights) # sum in new norm 
    b = np.sum(weights) # original integral over space treats each weight the same
    a = np.sum(A@LB@weights) # new integral and of smoothed (want a to equal b) 
    r = b/a
    weightsNew = r*LB@weights
    print("Weights Multiplied by " + str(r))
    featsVTKNames.append('TOTAL_AREA_OF_TISSUE')
    print("before " + str(np.sum(weights)))
    print("After normalization " + str(np.sum(r*A@LB@weights)))
    
    featsVTK.append(r*LB@weights)
    '''
    # add Tau in region per area 
    b = np.sum(A@feats[0]) # if wanted to preserve tau density 
    a = np.sum(A@LB@feats[0])
    r = b/a
    print("Tau density multiplied by " + str(r))
    featsNew = r*feats[0]
    featsVTK.append(r*LB@feats[0])
    print("before " + str(np.sum(A@feats[0])))
    print("After normalization " + str(np.sum(r*A@LB@feats[0])))
       
    featsVTKNames.append('TOTAL_TAU_PER_MMSQ')
    '''
    
    # add Tau total and normalize by that 
    totalTau = np.squeeze(feats[0])*np.squeeze(weights)
    b = np.sum(totalTau) # sum of total tau over space in varifold model
    a = np.sum(A@LB@totalTau)
    r = b/a
    print("Total tau multiplied by " + str(r))
    featsVTK.append(r*LB@totalTau)
    featsVTKNames.append('TOTAL_TAU')
    
    # add Tau total and divide by smoothed weights
    featsVTK.append((r*LB@totalTau)/weightsNew)
    featsVTKNames.append('TOTAL_TAU_PER_SMOOTHED_MM')
    
    vt.writeVarifoldVTK(YXZ,featsVTK,featsVTKNames,savename,polyData=polys)
    return

def writeVarifoldOntoMRI(varFile,pointMatFile,hdr,img,f0,f1,savename,down=1):
    params = sp.io.loadmat(pointMatFile)
    inds = np.asarray(params['x0x1x2'])
    I = ds.loadSeg(hdr,img)
    Im=np.asanyarray(I.dataobj)
    Idown=Im[0:Im.shape[0]:down,0:Im.shape[1]:down,0:Im.shape[2]:down]
    affineScale = np.eye(I.affine.shape[0])*down
    affineScale[-1,-1] = 1
    
    Inew = np.zeros_like(Idown)-1 # -1 if data not found 
    
    loc,weights,feats = loadVarifold(varFile) # assume in same order as in pointMatFile
    
    for i in range(weights.shape[0]):
        x0 = int(inds[i,0])
        x1 = int(inds[i,1])
        x2 = int(inds[i,2])
        
        Inew[x0,x1,x2] = feats[f0][i,f1]
    
    fileMap = nib.AnalyzeImage.make_file_map()
    fileMap['image'].fileobj = savename + '.img'
    fileMap['header'].fileobj = savename + '.hdr'
    totImage = nib.AnalyzeImage(Inew,I.affine*affineScale,file_map=fileMap)
    nib.save(totImage,savename + '.img')
    return