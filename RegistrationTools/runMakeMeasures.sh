#! /bin/bash

# Script that translates raw histology data (transformed) into discrete measure representation.

# Example: ./runMakeMeasures.sh -b 2 -o 1 -d 020922 -t 'Tau' -f 'entropy_round2' -x '_pil' -r 'redone_11222' -s 1

mriSuffix=''
mriFolder='itk'

while getopts b:o:d:t:f:x:r:s: flag
do
    case "${flag}" in
        b) brainNum=${OPTARG};; # 1 based 
        o) blockNum=${OPTARG};; # 1 based; total number of blocks 
        d) dateNew=${OPTARG};; # date for saving output
        t) stain=${OPTARG};; # stain for choosing Amyloid (0) or Tau (2)
        f) mriFolder=${OPTARG};; # entropy or itk or entropy_round2
        x) mriSuffix=${OPTARG};; # _pil
        r) redone=${OPTARG};; # redone_11222 for brain 2
        s) step=${OPTARG};;
    esac
done

echo $brainNum
echo $blockNum
echo $dateNew
echo $stain
echo $mriFolder
echo $mriSuffix
echo $step

case $stain in
    'Tau')
        stainNum=2
        ;;
    'Amyloid')
        stainNum=0
        ;;
esac

case $brainNum in
    2)
        numLab=15
        factor=4.4
        weightBy="2_2only_${mriSuffix}_${dateNew}"
        totSuff="_${mriSuffix}_${dateNew}"
        slicesToVTK="_FineDown4_jacobW${weightBy}"
        ;;
    5)
        numLab=16
        factor=1.7
        weightBy="2_5only_${mriSuffix}_${dateNew}"
        totSuff="_${mriSuffix}_${dateNew}"
        slicesToVTK="_FineDown4_jacobW${weightBy}"
        ;;   
esac

# Step 1: Estimate 2D Jacobian Transformation 
# Step 2: Downsample and Convert to Varifold + Select MTL Regions of Interest
# Step 3: Compute sliding window
# Step 4: write slices as VTK for visualization
# Step 5: Project to Surface (and Laplace Beltrami)
# Step 6: Make combined bar graph for 2 brains (regional densities)
# Step 7: Resample Downsampled via Gaussian 
case $step in
    1)
        python3 -c "import sys; sys.path.append('../'); sys.path.append('../VarifoldTools/'); import measureFunctions as mf; mf.getAreaWeightsExvivo($brainNum,[$blockNum-1],$stainNum,'$mriSuffix','$dateNew')"
        ;;
    2)
        python3 -c "import sys; sys.path.append('../'); sys.path.append('../VarifoldTools/'); import measureFunctions as mf; mf.downsampleExVivoAndConvert($brainNum,[0,1,2],$numLab,$stainNum,4,'$mriSuffix','$dateNew','$redone');mf.selectMTL($brainNum,'$stain','$mriSuffix','$dateNew')"
        python3 -c "import sys; sys.path.append('../'); sys.path.append('../VarifoldTools/'); import measureFunctions as mf;mf.selectMTL($brainNum,'$stain','$mriSuffix','$dateNew')"
        ;;
    3)
        python3 -c "import sys; sys.path.append('../'); sys.path.append('../VarifoldTools/'); import varifoldFunctions as vf; vf.makeBarFromVarifold($brainNum,'$stain',4,'$weightBy',lab2D=False,normalize=False,factor=1.0/$factor,ylim=30)"
        /cis/home/kstouff4/Documents/VarifoldTools/getSlidingWindow.sh $brainNum $totSuff
        ;;
    4) 
        cd /cis/home/kstouff4/Documents/SurfaceTools/
        ./writeSlicesAsVTK.sh $brainNum $blockNum $slicesToVTK
        ;;
    5)
        cd /cis/home/kstouff4/Documents/SurfaceTools/
        case $brainNum in
            2)
                amy=$redone/Surfaces/amygdala_toMai
                erc=$redone/Surfaces/ERC+Ext_label01_toMai_0803
                ca1=$redone/Surfaces/ca1_toMai
                sub=$redone/Surfaces/subiculum_toMai
                ;;
            5)
                amy=Surfaces/amygdala_toMai_1127
                erc=Surfaces/extension+erc_toMai_1127
                ca1=Surfaces/ca1_toMai_1127
                sub=Surfaces/subiculum_toMai_1127
                ;;
        esac
        matlab -nodisplay -r "load('/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain${brainNum}/3DSegmentations/${amy}.mat');sn='nothing';smootheLB(polys(:,2:end),YXZ,size(YXZ,1),0,sn,'/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain${brainNum}/3DSegmentations/${amy}_LBAll.mat');pause(10);quit"
        matlab -nodisplay -r "load('/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain${brainNum}/3DSegmentations/${erc}.mat');sn='nothing';smootheLB(polys(:,2:end),YXZ,size(YXZ,1),0,sn,'/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain${brainNum}/3DSegmentations/${erc}_LBAll.mat');pause(10);quit"
        matlab -nodisplay -r "load('/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain${brainNum}/3DSegmentations/${ca1}.mat');sn='nothing';smootheLB(polys(:,2:end),YXZ,size(YXZ,1),0,sn,'/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain${brainNum}/3DSegmentations/${ca1}_LBAll.mat');pause(10);quit"
        matlab -nodisplay -r "load('/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain${brainNum}/3DSegmentations/${sub}.mat');sn='nothing';smootheLB(polys(:,2:end),YXZ,size(YXZ,1),0,sn,'/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain${brainNum}/3DSegmentations/${sub}_LBAll.mat');pause(10);quit"

        cd /cis/home/kstouff4/Documents/RegistrationTools
        python3 -c "import sys; sys.path.append('../'); sys.path.append('../VarifoldTools/'); import measureFunctions as mf; mf.sampleToSurfaces($brainNum,'$stain','$mriSuffix','$dateNew','$redone'); mf.writeSurfacesAndLB($brainNum,'$stain','$mriSuffix','$dateNew','$redone')"
        ;;
    6)
        python3 -c "import sys; sys.path.append('../'); sys.path.append('../VarifoldTools/'); import varifoldFunctions as vf; vf.makeDoubleBarFromVarifold([2,5],'$stain',4,['2_2only_${mriSuffix}_${dateNew}','2_5only_${mriSuffix}_${dateNew}'],lab2D=False,normalize=False,factors=[1.0/4.4,1.0/1.7],ylim=30)"
        ;;
    7)
        python3 -c "import sys; sys.path.append('../'); sys.path.append('../VarifoldTools/'); import measureFunctions as mf;mf.selectFour($brainNum,'$stain','$mriSuffix','$dateNew');mf.resampleIsotropicGauss($brainNum,'$stain','$mriSuffix','$dateNew')"
        ;;
esac

 

