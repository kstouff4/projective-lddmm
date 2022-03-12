#! /bin/bash

# Script that outlines all steps
# Example: ./runFullPipeline.sh -b 2 -o 1 -s 0 -d 020222 -i 500 -t 'Tau' -f 'entropy' -x '_pil'

mriSuffix=''
mriFolder='itk'

while getopts b:o:s:i:d:t:f:x: flag
do
    case "${flag}" in
        b) brainNum=${OPTARG};; # 1 based 
        o) blockNum=${OPTARG};; # 1 based 
        s) step=${OPTARG};; # 0 = estimating ST and 3D transformations, 1 = estimating 3D and 2D 
        i) iter=${OPTARG};; # number of iterations to repeat steps 3+4
        d) dateNew=${OPTARG};; # date for saving output
        t) stain=${OPTARG};; # stain for choosing Amyloid (0) or Tau (2)
        f) mriFolder=${OPTARG};; # entropy or itk or entropy_round2
        x) mriSuffix=${OPTARG};; # _pil
    esac
done

echo $brainNum
echo $blockNum
echo $step
echo $iter
echo $dateNew
echo $stain
echo $mriFolder
echo $mriSuffix

# Step 0: Estimate 3D Transformations without scattering transform (1000 iter)
# Args: BrainNum, blockNum
cd /cis/home/kstouff4/Documents/HistoMRIPipeline
case $step in 

    0)
        matlab -nodisplay -r "getMRICoords_by_katie('Brain$brainNum','$blockNum');align_v01_diffeo_by_katie_single('Brain$brainNum','$blockNum','$stain',$iter,false);quit()"
        
        cd /cis/home/kstouff4/Documents/RegistrationTools
        
        mkdir /cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain$brainNum/histology/alignment/block$blockNum/$mriFolder/
        python -c "import biasCorrect as bc; fname='/cis/home/kstouff4/Documents/HistoMRIPipeline/subject_Brain'+str($brainNum) + '_block_' + str($blockNum) + '_alignment_diffeo_test9_by_katie_single/' + str($iter) + '_iters/mris.mat'; bc.runBias($brainNum,$blockNum,'$mriFolder','$mriSuffix',fname)"

        # Step 1: Compute and Cache Scattering Transform of Histology Images
        # Args: BrainNum, blockNum (fixed paths)

        python3 -c "import ScatteringFunctions as sf; sf.preProcessScatter($brainNum, $blockNum, '$stain')"

        # Step 2: Estimate f(): ST --> grayscale by also estimating 2D Transformations (500 iter); fix f()
        mkdir /cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain$brainNum/histology/down_000/AD_Hip$blockNum/outputFigs/F_$dateNew/
        mkdir /cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain$brainNum/histology/down_000/AD_Hip$blockNum/outputFigs/F_$dateNew/$stain
        python3 -c "import Estimate2DTransforms as e2t; paramsf = '/cis/home/kstouff4/Documents/HistoMRIPipeline/subject_Brain'+str($brainNum) + '_block_' + str($blockNum) + '_alignment_diffeo_test9_by_katie_single/' + str($iter) + '_iters/relTrans.mat'; dMRI = e2t.getResolution(paramsf);e2t.estimateF($brainNum,$blockNum,'$stain','$mriFolder','$mriSuffix','F_$dateNew',dMRI)"
        
        ;;
    1)

        # 3D and 2D Transformations are re-Initialized to Id 
        # Step 3: Estimate 3D Transformations with fixed scattering transform (500 iter)
        # Args: BrainNum, blockNum (either fAI or fI)
        
        ##matlab -nodisplay -r "align_v01_diffeo_by_katie_single_round2('Brain$brainNum','$blockNum','$stain',500,'fAI_F_$dateNew.tif','$dateNew');quit()"

        cd /cis/home/kstouff4/Documents/RegistrationTools
        ##mkdir /cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain$brainNum/histology/alignment/block$blockNum/${mriFolder}_round2/
        ##python -c "import biasCorrect as bc; fname='/cis/home/kstouff4/Documents/HistoMRIPipeline/subject_Brain'+str($brainNum) + '_block_' + str($blockNum) + '_alignment_diffeo_test9_by_katie_single_round2/500_iters_$dateNew/mris2.mat'; bc.runBias($brainNum,$blockNum,'$mriFolder'+'_round2','$mriSuffix',fname)"
        
        # Step 4: Estimate 2D Transformations with fixed 3D Transformations (200 iter)
        # Args: BrainNum, blockNum
        mkdir /cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain$brainNum/histology/down_000/AD_Hip$blockNum/outputFigs/$dateNew/
        mkdir /cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain$brainNum/histology/down_000/AD_Hip$blockNum/outputFigs/$dateNew/$stain
        python3 -c "import Estimate2DTransforms as e2t; paramsf = '/cis/home/kstouff4/Documents/HistoMRIPipeline/subject_Brain'+str($brainNum) + '_block_' + str($blockNum) + '_alignment_diffeo_test9_by_katie_single_round2/500_iters_$dateNew/relTrans2.mat'; dMRI = e2t.getResolution(paramsf);  e2t.estimate2DPhi($brainNum,$blockNum,'$stain','$mriFolder'+'_round2','$mriSuffix','$dateNew','F_$dateNew',dMRI)"


        # Step 5: Return to Step 3 with current 3D and update as needed 
        ;;
esac