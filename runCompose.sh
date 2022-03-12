#! /bin/bash

# Script that composes transformations following runFullPipeline.sh

# Example Full Block: ./runCompose.sh -b 2 -o 1 -d 020922 -t 'Tau' -f 'entropy_round2' -x '_pil'
# Example Single Slice: ./runCompose.sh -b 2 -o 2 -d 021222 -t 'Tau' -f 'entropy_round2' -x '_pil' -m '/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/histology/down_000/AD_Hip2/outputFigs/021222/Tau/Brain2_Block2_L5.mat' -i 0

mriSuffix=''
mriFolder='itk'
suffixSave=''

while getopts b:o:d:t:f:x:m:i:s: flag
do
    case "${flag}" in
        b) brainNum=${OPTARG};; # 1 based 
        o) blockNum=${OPTARG};; # 1 based 
        d) dateNew=${OPTARG};; # date for saving output
        t) stain=${OPTARG};; # stain for choosing Amyloid (0) or Tau (2)
        f) mriFolder=${OPTARG};; # entropy or itk or entropy_round2
        x) mriSuffix=${OPTARG};; # _pil
        m) matFile=${OPTARG};; # matfile if not none for manual alignment
        i) sliceNum=${OPTARG};; # slice index (0 based) for manual alignment 
        s) suffixSave=${OPTARG};; # saving suffix (e.g. XM if manual)
    esac
done

echo $brainNum
echo $blockNum
echo $dateNew
echo $stain
echo $mriFolder
echo $mriSuffix
echo $matFile
echo $sliceNum
echo $suffixSave

if [ -z "$sliceNum" ]
then
      python3 -c "import FullPipelineFunctions as fpf; fpf.getComposition($brainNum,$blockNum,'$stain','$dateNew','$mriSuffix')"
      case $brainNum in 
        2)
            ./runGetInfo.sh $brainNum $blockNum "_${mriSuffix}_${dateNew}.npz" redone_11222 "_${mriSuffix}_${dateNew}"
            ;;
        5)
            ./runGetInfo.sh $brainNum $blockNum "_${mriSuffix}_${dateNew}.npz" redone_1105 "_${mriSuffix}_${dateNew}"
            ;;
      esac 
else  
      if [ -z "$suffixSave" ]
      then
          suffixSave=''
      fi
      python3 -c "import FullPipelineFunctions as fpf; fpf.getCompositionManual($brainNum,$blockNum,$sliceNum,'$stain','$dateNew','$mriSuffix','$matFile')"
      if [[ "$matFile" != *"XM"* ]]; then
          echo "getting labels"
          case $brainNum in
              2)
                  locs=(5 6 7 8 9 10 11 12 13 14 15 16 17 18 19)
                  python3 -c "import FullPipelineFunctions as fpf; fpf.getInfo($brainNum, 2, $blockNum-1, $sliceNum, ${locs[$sliceNum]},suffix='_${mriSuffix}_${dateNew}.npz',redone='redone_11222',suffixSave='_${mriSuffix}_${dateNew}${suffixSave}')"
                  ;;
              5)
                  locs=(1 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75)
                  if (($blockNum == 3))
                  then
                      locs=(3 7 12 17 22 27 32 37)
                  fi
                  python3 -c "import FullPipelineFunctions as fpf; fpf.getInfo($brainNum, 2, $blockNum-1, $sliceNum, ${locs[$sliceNum]},suffix='_${mriSuffix}_${dateNew}.npz',redone='redone_1105',suffixSave='_${mriSuffix}_${dateNew}${suffixSave}')"
                  ;;
          esac
      fi
fi

#for (( i=2; i<=$blockNum; i++ ))
#do  
#   python3 -c "import FullPipelineFunctions as fpf; fpf.getComposition($brainNum,$i,'$stain','$dateNew','$mriSuffix')"
#done

# Get labels based on which folder has most up to date labels per brain 

    
    