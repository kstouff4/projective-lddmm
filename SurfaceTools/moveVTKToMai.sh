#! /bin/bash

python3 - << EOF

import vtkFunctions as vt
import glob

files1 = glob.glob('/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain5/Block1/Tau/3DVarifold/VTKPartFiles/*5only.vtk')
files2 = glob.glob('/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain5/Block2/Tau/3DVarifold/VTKPartFiles/*5only.vtk')
files3 = glob.glob('/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain5/Block3/Tau/3DVarifold/VTKPartFiles/*5only.vtk')
matfile = '/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain5/Params/A_toMai_1127.mat'

for f in files1:
    sn = f.split('.')[0] + '_inMRI.vtk'
    vt.moveToMai(matfile,f,sn,inverse=1,shift=[30.0,25.0,17.5])

for f in files2:
    sn = f.split('.')[0] + '_inMRI.vtk'
    vt.moveToMai(matfile,f,sn,inverse=1,shift=[30.0,25.0,17.5])

for f in files3:
    sn = f.split('.')[0] + '_inMRI.vtk'
    vt.moveToMai(matfile,f,sn,inverse=1,shift=[30.0,25.0,17.5])

EOF