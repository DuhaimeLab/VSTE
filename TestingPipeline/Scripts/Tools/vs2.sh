#!/bin/bash

echo start

rm -r "VS2_${1}"
mkdir "VS2_${1}"
cd "VS2_${1}"
pwd

samplename="/scratch/duhaimem_root/duhaimem/shared_data/VSTE/TestingSet/metagenomic_testing_set_${1}.fna"
threads=$2
path="/scratch/duhaimem_root/duhaimem/shared_data/VSTE/RunningToolsOnTestingSet/VS2/VS2_${1}"
threads=10

virsorter run -w $path -i $samplename --include-groups all -j $threads all

echo done
