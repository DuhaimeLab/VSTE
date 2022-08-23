#!/bin/bash

echo start
echo $1
samplename="/scratch/duhaimem_root/duhaimem/shared_data/VSTE/TestingSet/metagenomic_testing_set_${1}.fna"
database="/home/hegartyb/checkv-db-v0.6/"
threads=$2
echo $threads

rm -r "Checkv_${1}"
mkdir "Checkv_${1}"
cd "Checkv_${1}"

echo contamination
checkv contamination $samplename ./ -t $threads -d $database
echo completeness
checkv completeness $samplename ./ -t $threads -d $database
echo complete genomes
checkv complete_genomes $samplename ./
echo quality summary
checkv quality_summary $samplename ./

echo done
