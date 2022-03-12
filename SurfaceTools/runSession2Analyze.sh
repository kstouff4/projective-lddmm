#! /bin/bash
# Script to get individual image segmentations out of seg3D session, save them, get a mask, and get individual surfaces from all the labeled tissue and the ERC separately. 
# Example: ./runSession2Analyze.sh /cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain4/claire_update.ses /cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain4/3DSegmentations/redone_2822/ '\(all\|AD_Hip\)' 3

# smoothed imagge only used for interpolation downstream of deforming labels onto 2D histology
cd /cis/home/kstouff4/Documents/SurfaceTools/
rm -rf /tmp/session2*/

matlab -nodisplay -r "session2Analyze('$1');pause(10);quit"

mkdir $2 # ensure directory is made 
cp *.img $2
cp *.hdr $2
rm *.img
rm *.hdr 

#matlab -nodisplay -r "session2Analyze('/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain3/3DSegmentations/post-Norway_3.ses');pause(10);quit"

#matlab -nodisplay -r "session2Analyze('/cis/home/kstouff4/Documents/HistoMRIPipeline/Brain2/Block1/Tau/eXuSegmentations/seg_updated_5.ses');pause(10);quit"

#cp *.img /cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/2DSegmentations/
#cp *.hdr /cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/2DSegmentations/
#rm *.img 
#rm *.hdr

cd $2
find . -type f -name "* *" | while read file; do mv "$file" ${file// /_}; done

# get list of files
y2=$(find . -maxdepth 1 -iregex ".*\(erc\|tec\).*img")
if [ ! "$y2" ];then
    OUT2=''
fi

for I in $y2
do
    OUT2=${OUT2:+$OUT2,}\'$2${I:2}\'
done
echo $OUT2

y=$(find . -maxdepth 1 -name "*img*" ! -iregex ".*$3.*")
for I in $y
do
    OUT=${OUT:+$OUT,}\'$2${I:2}\'
done

echo $OUT

cd /cis/home/kstouff4/Documents/HistoMRIPipeline/seg2Surf/

case $4 in
    2)
        #python3 -c "import glob; import numpy as np; import segTools as st; f = [$OUT]; f.sort(key=lambda x: x.lower()); st.getAll2(f,'$2all.img'); f2 = [$OUT2]; st.getAll2(f2,'$2TEC+ERC.img'); st.allToMask('$2all.img'); st.mergeSeg2(f,'$2allMerge')"
        #;;
        python3 -c "import glob; import numpy as np; import segTools as st; f = [$OUT]; f.sort(key=lambda x: x.lower()); st.getAll2(f,'$2all.img'); st.allToMask('$2all.img'); st.mergeSeg2(f,'$2allMerge')"
        ;;
        
    3)
        python3 -c "import glob; import numpy as np; import segTools as st; f = [$OUT]; f.sort(key=lambda x: x.lower()); print('files are:'); print(f); st.getAll(f,'$2all.img'); st.allToMask('$2all.img'); st.mergeSeg(f,'$2allMerge');st.smootheWithMode('$2allMerge.hdr','$2allMerge.img',1,'$2allMerge_smoothe.hdr','$2allMerge_smoothe.img')"
        ;;
esac 

#/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain3/3DSegmentations/redone_1105/all.img");st.allToMask("/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain3/3DSegmentations/redone_1105/all.img");st.mergeSeg(f,"/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain3/3DSegmentations/redone_1105/allMerge"); st.smootheWithMode("/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain3/3DSegmentations/redone_1105/allMerge.hdr","/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain3/3DSegmentations/redone_1105/allMerge.img",1,"/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain3/3DSegmentations/redone_1105/allMerge_smoothe.hdr","/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain3/3DSegmentations/redone_1105/allMerge_smoothe.img")'

# Combine ERC + Extension
#python3 -c 'import glob; import numpy as np; import segTools as st; f = ["/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain5/3DSegmentations/redone_1105/erc Label 1.img", "/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain5/3DSegmentations/redone_1105/claire extension Label 1.img"];st.getAll(f,"/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain5/3DSegmentations/redone_1105/extension+erc.img")'

#cd /cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/3DSegmentations/Surfaces

#/cis/home/kstouff4/Documents/HistoMRIPipeline/seg2Surf/seg2Surf.sh "/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain7/3DSegmentations/all_mask.img" 0.75 0.125 1

#/cis/home/kstouff4/Documents/HistoMRIPipeline/seg2Surf/seg2Surf.sh "/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain7/3DSegmentations/erc_Label_1.img" 0.75 0.125 1
