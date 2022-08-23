#!/bin/bash

echo start

rm "Vibrant_${1}"
mkdir "Vibrant_${1}"
cd "Vibrant_${1}"
pwd

samplename="/scratch/duhaimem_root/duhaimem/shared_data/VSTE/TestingSet/metagenomic_testing_set_${1}.fna"
threads=$2
database="/home/hegartyb/VIBRANT/databases/"
vibrantloc="/home/hegartyb/VIBRANT/files/"

python3 /home/hegartyb/VIBRANT/VIBRANT_run.py -i $samplename -t $threads -folder ./ -d $database -m $vibrantloc

echo done
