#!/bin/bash

echo start

rm -r "DVF_{1}"
mkdir "DVF_${1}"

samplename="/scratch/duhaimem_root/duhaimem/shared_data/VSTE/TestingSet/metagenomic_testing_set_${1}.fna"
threads=$2
outfile="/scratch/duhaimem_root/duhaimem/shared_data/VSTE/RunningToolsOnTestingSet/DVF/DVF_${1}"

python /scratch/kwigg_root/kwigg/hegartyb/DeepVirFinder/dvf.py -i $samplename -o $outfile -c $threads

echo done
