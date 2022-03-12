#! /bin/bash
# Script to get smoothed version of functions using LB

# Brain 2 Params:
# sumCorticalColumns_Template_0803_Amyloid.mat
# erc_Template_0803.mat
# amygdala = 1, erc = 5, subiculum = 11

# Brain 5 Params:
# sumCorticalColumns_Template_0803_Tau.mat
# erc_Template_0803.mat

# run LB Basis
matlab -nodisplay -r "load('/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain3/3DSegmentations/Surfaces/amygdala_toMai_0803.mat');sn='nothing';disp(size(YXZ,1));smootheLB(polys(:,2:end),YXZ,size(YXZ,1),0,sn,'/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain3/3DSegmentations/Surfaces/amygdala_toMai_0803_LBAll.mat');pause(10);quit"

# alter Varifolds
cd /cis/home/kstouff4/Documents/VarifoldTools/

# run with -1 as third argument if just 1 region
python3 -c 'import varifoldFunctions as var; import scipy as sp; import numpy as np; mf = "/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain3/3DSegmentations/Surfaces/amygdala_toMai_0803.mat"; bf = "/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain3/3DSegmentations/Surfaces/amygdala_toMai_0803_LBAll.mat"; params2 = sp.io.loadmat(bf); B = np.asarray(params2["B"]); D = np.asarray(params2["D"]); A = np.asarray(params2["A"]); varFile = "/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain3/MultiBlock/Tau/3DVarifold/VarifoldPartFiles/Amygdala_FineDown4_NonejacobW2_NN_Surface_down4_3only.npz"; sn = "/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain3/MultiBlock/Tau/3DVarifold/VarifoldVTKFiles/Amygdala_FineDown4_NonejacobW2_NN_Surface_down4_3only_LBAll_k2.vtk"; var.writeRegionVarifoldToVTK(varFile,mf,-1,"ERC",sn,B,D,A,2)'

# run with -1 as third argument if just 1 region
#python3 -c 'import varifoldFunctions as var; import scipy as sp; import numpy as np; mf = "/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain5/3DSegmentations/Surfaces/extension+erc_toMai_0803.mat"; varFile = "/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain2/MultiBlock/Tau/3DVarifold/VarifoldPartFiles/Extension+ERC_FineDown4_NonejacobW2_NN_Surface_down4_2only.npz"; sn = "/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain2/MultiBlock/Tau/3DVarifold/VarifoldVTKFiles/Extension+ERC_FineDown4_NonejacobW2_NN_Surface_down4_2only.vtk"; var.writeRegionVarifoldToVTK(varFile,mf,0,"ERC",sn,None,None,None,0)'


# run heat Eq
#matlab -nodisplay -r "load('/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain2/Thickness/sumCorticalColumns_Template_0803_Top_0912_Tau.mat');load('/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain2/Thickness_flipped/erc_Template_0803.mat');iter=25;eps=0.1;heatEq(faces+1,vertices,iter,eps,[feats0 weights],'/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain2/Thickness/sumCorticalColumns_Template_0803_Top_0912_Tau_HQ25.mat');pause(10);quit"

