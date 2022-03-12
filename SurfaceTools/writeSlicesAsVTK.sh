#! /bin/bash

brainNum=$1
blockNum=$2
baseSuff=$3

python3 - << EOF
import vtkFunctions as vtkF
import sys
sys.path.append('../')
sys.path.append('../VarifoldTools/')
import varifoldFunctions as vF
import numpy as np

baseSuffix = "$baseSuff.npz" #"_FineDown4_jacobW2_3only.npz"
if ($brainNum == 2):
    baseIn = ["L5", "L6", "L7","L8", "L9", "L10", "L11", "L12", "L13", "L14", "L15", "L16", "L17", "L18", "L19"]
    if ($blockNum == 2):
        baseIn = ["L5","L6", "L7","L8", "L9", "L10", "L11", "L12", "L13", "L14", "L15", "L16", "L17", "L18", "L19"]
    if ($blockNum == 3):
        baseIn = ["L5", "L6", "L7","L8", "L9", "L10", "L11"] #"L12", "L13", "L14", "L15", "L16", "L17", "L18", "L19"]
elif ($brainNum == 5):
    if ($blockNum == 1):
        baseIn = ["L1","L5","L10","L15","L20","L25","L30","L35","L40","L45","L50","L55","L60"]
    elif ($blockNum == 2):
        baseIn = ["L1","L5","L10","L15","L20","L25","L30","L35","L40","L45","L50","L55","L65","L70"]
    elif ($blockNum == 3):
        baseIn = ["L3", "L7","L12","L17","L22","L27","L32","L37"]
baseName = "/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain${brainNum}/Block${blockNum}/Tau/3DVarifold/VarifoldPartFiles/"
vtkBaseName = "/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain${brainNum}/Block${blockNum}/Tau/3DVarifold/VTKPartFiles/"
vtkBaseSuffix = "$baseSuff.vtk"

for b in baseIn:
    loc,weight,feats = vF.loadVarifold(baseName + b + baseSuffix)
    # Assume the loc is in r x c x 1 x 3
    loc = np.squeeze(loc)
    weight = np.squeeze(weight)
    print(loc.shape)
    print(weight.shape)
    r = loc.shape[0]
    c = loc.shape[1]
    locNew = np.zeros((2*r + 2*c,3))
    locNew[0:r] = loc[:,0,...]
    locNew[r:r+c] = loc[-1,:,...]
    locNew[r+c:r+r+c] = loc[:,-1,...]
    locNew[r+r+c:] = loc[0,:,...]
    savename = vtkBaseName + b + vtkBaseSuffix
    weightsNew = np.zeros((2*r + 2*c,1))
    weightsNew[0:r] = weight[:,0,None]
    weightsNew[r:r+c] = weight[-1,:,None]
    weightsNew[r+c:r+r+c] = weight[:,-1,None]
    weightsNew[r+r+c:] = weight[0,:,None]
    vtkF.writeVTK(locNew,[weightsNew],['AREA'],savename,polyData=None)

EOF
