#! /bin/bash

suff=$2 #'__pil_020922'
brainNum=$1 # 2

python3 - << EOF

import sys
sys.path.append('../')
sys.path.append('../VarifoldTools/')
import varifoldFunctions as vf
import numpy as np
import matplotlib.pyplot as plt

width = 1
center = '2'
if ($brainNum == 5):
    sub_5 = "/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain5/MultiBlock/Tau/3DVarifold/VarifoldPartFiles/Subiculum_FineDown4_jacobW2_5only$suff.npz"
    ca1_5 = "/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain5/MultiBlock/Tau/3DVarifold/VarifoldPartFiles/CA1_FineDown4_jacobW2_5only$suff.npz"
    amy_5 = "/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain5/MultiBlock/Tau/3DVarifold/VarifoldPartFiles/Amygdala_FineDown4_jacobW2_5only$suff.npz"
    erc_5 = "/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain5/MultiBlock/Tau/3DVarifold/VarifoldPartFiles/ERC+Ext_FineDown4_jacobW2_5only$suff.npz"
    all_5 = "/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain5/MultiBlock/Tau/3DVarifold/VarifoldPartFiles/Hippocampus_Sub_FineDown4_jacobW2_5only$suff.npz"
    savebase5 = "/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain5/MultiBlock/Tau/BarAndLine/vertPermm$suff"
    inds1, ca1, ca1A = vf.plotSlidingWindow([ca1_5],width,savebase5+'_ca1_' + str(width) + '.png',2,1,-1,45,0)
    inds1, sub, subA = vf.plotSlidingWindow([sub_5],width,savebase5+'_sub_' + str(width) + '.png',2,1,-1,45,0)
    inds1, amy, amyA = vf.plotSlidingWindow([amy_5],width,savebase5+'_amy_' + str(width) + '.png',2,1,-1,45,0)
    inds1, erc, ercA = vf.plotSlidingWindow([erc_5],width,savebase5+'_erc_' + str(width) + '.png',2,1,-1,45,0)
    inds1, hippo, hippoA = vf.plotSlidingWindow([all_5],width,savebase5+'_hippo_' + str(width) + '.png',2,1,-1,45,0)
    vf.plotTrendsSeaborn(savebase5+'bar+trend_amyErc.png',inds1,[amy,erc],[amyA,ercA],legLab=["Amygdala","ERC"],title="Tau Tangle Density Rostral to Caudal",xlab="Center of " + center + " mm window (mm, Mai Z Axis)",ylab="Tau Tangles / mm$^2$",cols=['green','plum'],trend=False,factor=1/1.7,xlim=30)
    vf.plotTrendsSeaborn(savebase5+'bar+trend_ca1Sub.png',inds1,[ca1,sub],[ca1A,subA],legLab=["CA1","Subiculum"],title="Tau Tangle Density Rostral to Caudal",xlab="Center of " + center + " mm window (mm, Mai Z Axis)",ylab="Tau Tangles / mm$^2$",cols=['lime','fuchsia'],trend=True,factor=1/1.7,xlim=30)

    vf.plotTrendsSeaborn(savebase5+'bar+trend_hippo_Sub.png',inds1,[hippo],[hippoA],legLab=["Hippocampus"],title="Tau Tangle Density Rostral to Caudal",xlab="Center of " + center + " mm window (mm, Mai Z Axis)",ylab="Tau Tangles / mm$^2$",cols=['blue'],trend=True,factor=1/1.7,xlim=30)
    vf.plotTrendsSeaborn(savebase5+'bar_hippo_Sub.png',inds1,[hippo],[hippoA],legLab=["Hippocampus"],title="Tau Tangle Density Rostral to Caudal",xlab="Center of " + center + " mm window (mm, Mai Z Axis)",ylab="Tau Tangles / mm$^2$",cols=['blue'],trend=False,factor=1/1.7,xlim=30)
    vf.plotTrendsSeaborn(savebase5+'bar_hippo_Sub+Amygdala.png',inds1,[amy,hippo],[amyA,hippoA],legLab=["Amygdala","Hippocampus"],title="Tau Tangle Density Rostral to Caudal",xlab="Center of " + center + " mm window (mm, Mai Z Axis)",ylab="Tau Tangles / mm$^2$",cols=['green','blue'],trend=False,factor=1/1.7,xlim=30)
    

elif ($brainNum == 2):
    #_,_,_,_ = vf.plotMaiCoronalTau([ca1_5,sub_5,amy_5,erc_5],width,savebase5+'_4struct_' + str(width) + '_slice5.png',2,5,0)
    sub_2 = "/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain2/MultiBlock/Tau/3DVarifold/VarifoldPartFiles/Subiculum_FineDown4_jacobW2_2only$suff.npz"
    ca1_2 = "/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain2/MultiBlock/Tau/3DVarifold/VarifoldPartFiles/CA1_FineDown4_jacobW2_2only$suff.npz"
    amy_2 = "/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain2/MultiBlock/Tau/3DVarifold/VarifoldPartFiles/Amygdala_FineDown4_jacobW2_2only$suff.npz"
    erc_2 = "/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain2/MultiBlock/Tau/3DVarifold/VarifoldPartFiles/ERC+Ext_FineDown4_jacobW2_2only$suff.npz"
    savebase2 = "/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain2/MultiBlock/Tau/BarAndLine/vertPermm$suff"
    all_2 = "/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain2/MultiBlock/Tau/3DVarifold/VarifoldPartFiles/Hippocampus_Sub_FineDown4_jacobW2_2only$suff.npz"
    inds2, ca1, ca1A = vf.plotSlidingWindow([ca1_2],width,savebase2+'_ca1_' + str(width) + '.png',2,1,-5,50,0)
    inds2, sub, subA = vf.plotSlidingWindow([sub_2],width,savebase2+'_sub_' + str(width) + '.png',2,1,-5,50,0)
    inds2, amy, amyA = vf.plotSlidingWindow([amy_2],width,savebase2+'_amy_' + str(width) + '.png',2,1,-5,50,0)
    inds2, erc, ercA = vf.plotSlidingWindow([erc_2],width,savebase2+'_erc_' + str(width) + '.png',2,1,-5,50,0)
    inds2, hippo, hippoA = vf.plotSlidingWindow([all_2],width,savebase2+'_hippo_' + str(width) + '.png',2,1,-5,50,0)
    #_,_,_,_ = vf.plotMaiCoronalTau([ca1_2,sub_2,amy_2,erc_2],width,savebase2+'_4struct_' + str(width) + 'slice5.png',2,5,0)
    vf.plotTrendsSeaborn(savebase2+'bar+trend_amyErc.png',inds2,[amy,erc],[amyA,ercA],legLab=["Amygdala","ERC"],title="Tau Tangle Density Rostral to Caudal",xlab="Center of " + center + " mm window (mm, Mai Z Axis)",ylab="Tau Tangles / mm$^2$",cols=['green','plum'],trend=False,factor=1/4.4,xlim=30)
    vf.plotTrendsSeaborn(savebase2+'bar+trend_ca1Sub.png',inds2,[ca1,sub],[ca1A,subA],legLab=["CA1","Subiculum"],title="Tau Tangle Density Rostral to Caudal",xlab="Center of " + center + " mm window (mm, Mai Z Axis)",ylab="Tau Tangles / mm$^2$",cols=['lime','fuchsia'],trend=True,factor=1/4.4,xlim=30)

    vf.plotTrendsSeaborn(savebase2+'bar+trend_hippo_Sub.png',inds2,[hippo],[hippoA],legLab=["Hippocampus"],title="Tau Tangle Density Rostral to Caudal",xlab="Center of " + center + " mm window (mm, Mai Z Axis)",ylab="Tau Tangles / mm$^2$",cols=['blue'],trend=True,factor=1/4.4,xlim=30)
    vf.plotTrendsSeaborn(savebase2+'bar_hippo_Sub.png',inds2,[hippo],[hippoA],legLab=["Hippocampus"],title="Tau Tangle Density Rostral to Caudal",xlab="Center of " + center + " mm window (mm, Mai Z Axis)",ylab="Tau Tangles / mm$^2$",cols=['blue'],trend=False,factor=1/4.4,xlim=30)
    vf.plotTrendsSeaborn(savebase2+'bar_hippo_Sub+Amygdala.png',inds2,[amy,hippo],[amyA,hippoA],legLab=["Amygdala","Hippocampus"],title="Tau Tangle Density Rostral to Caudal",xlab="Center of " + center + " mm window (mm, Mai Z Axis)",ylab="Tau Tangles / mm$^2$",cols=['green','blue'],trend=False,factor=1/4.4,xlim=30)
    
    #f,ax = plt.subplots(figsize=(6,10))
#minVal = np.min(subCA1)+0.001
#maxVal = np.max(amyERC) 
#amyERCRat = amyERC/maxVal
#subCA1Rat = subCA1/maxVal
#ax.plot(amyERC,inds1,c="r",label="Amygdala+ERC")
#ax.plot(subCA1,inds1,c="b",label="Subiculum+CA1")
#plt.legend()
#ax.set_xlabel("Ratio of Tau Tangles in 2 mm window to Max along Z axis") 
#ax.set_ylabel("Center of 2 mm window") f.savefig("/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain5/MultiBlock/Tau/BarAndLine/amyERCsubCA1ratios_max2_vert_permm_raw.png")

#varFiles = ["/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain2/MultiBlock/Tau/3DVarifold/VarifoldPartFiles/Subiculum_FineDown4_NonejacobW2_NN_Surface_down4_2only.npz","/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain2/MultiBlock/Tau/3DVarifold/VarifoldPartFiles/CA1_FineDown4_NonejacobW2_NN_Surface_down4_2only.npz"]
#savename = "/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain2/MultiBlock/Tau/BarAndLine/sub+ca1_slidingWin2_vert_permm.png" inds, subCA1_2 = vf.plotSlidingWindow(varFiles,width,savename,2,1,-5,50,0) 
#varFiles = ["/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain2/MultiBlock/Tau/3DVarifold/VarifoldPartFiles/Amygdala_FineDown4_NonejacobW2_NN_Surface_down4_2only.npz","/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain2/MultiBlock/Tau/3DVarifold/VarifoldPartFiles/Extension+ERC_FineDown4_NonejacobW2_NN_Surface_down4_2only.npz"]
#savename = "/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain2/MultiBlock/Tau/BarAndLine/amy+ext+erc_slidingWin2_vert_permm.png"
#inds, amyERC_2 = vf.plotSlidingWindow(varFiles,width,savename,2,1,-5,50,0) 

#f,ax = plt.subplots(figsize=(10,6)) 
#f2,ax2 = plt.subplots(figsize=(10,6)) 

#minVal_2 = np.min(subCA1_2)+0.001 
#maxVal_2 = np.max(amyERC_2) 
#amyERCRat_2 = amyERC_2/maxVal_2 
#subCA1Rat_2 = subCA1_2/maxVal_2 

#ax2.plot(inds2,amy2,c="green",label="Amygdala Brain 2") 
#ax2.plot(inds2,erc2,c="purple",label="ERC Brain 2")
#ax2.plot(inds1,erc,c="plum",label="ERC Brain 5")
#ax.plot(inds2,sub2,c="fuchsia",label="Subiculum Brain 2")
#ax.plot(inds1,sub,c="pink",label="Subiculum Brain 5")
#ax2.plot(inds1,amy,c="olive",label="Amygdala Brain 5") 
#ax.plot(inds1,ca1,c="yellow",label="CA1 Brain 5")
#ax.plot(inds2,ca12,c="lime",label="CA1 Brain 2")
#ax.legend()
#ax.set_ylabel("Tau Tangle Density in 2 mm window along Z axis") 
#ax.set_xlabel("Center of 2 mm window") 
#ax2.set_ylabel("Tau Tangle Density in 2 mm window to Max along Z axis") 
#ax2.set_xlabel("Center of 2 mm window")
#ax2.legend()
#f.savefig("/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain5/MultiBlock/Tau/BarAndLine/subCA1_max2_2and5_permm_raw_" + str(width) + ".png")
#f2.savefig("/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain5/MultiBlock/Tau/BarAndLine/amyERC_max2_2and5_permm_raw_" + str(width) + ".png")
#f,ax = plt.subplots(figsize=(10,6))
#ax.plot(inds1,all5,c="blue",label="Brain 5")
#ax.plot(inds2,all2,c="red",label="Brain 2")
#ax.set_xlabel("Center of " + str(2*width) + " mm window")
#ax.set_ylabel("Tau Tangle Density per histological area")
#ax.legend()
#f.savefig("/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain5/MultiBlock/Tau/BarAndLine/4struct_2and5_permm_raw_" + str(width) + ".png")

EOF