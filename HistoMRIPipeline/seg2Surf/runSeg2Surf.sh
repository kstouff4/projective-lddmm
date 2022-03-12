#!/bin/bash

cd /cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/3DSegmentations/redone_11222/Surfaces
#cd /cis/home/kstouff4/Documents/datasets/ADNI_20-20-20/MCI/1117/
# 0 if don't shift (only use for Brain 2)

/cis/home/kstouff4/Documents/HistoMRIPipeline/seg2Surf/seg2Surf.sh "/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/3DSegmentations/redone_11222/amygdala_Label.img" 0.75 0.125 1

#/cis/home/kstouff4/Documents/HistoMRIPipeline/seg2Surf/seg2Surf.sh "/cis/home/skulason/Documents/ADNI_reorganized/data/mci/manseg/seg/016_S_1117_mo06_ERC_and_TEC.img" 0.75 1.0 1