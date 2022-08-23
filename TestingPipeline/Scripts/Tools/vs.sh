#!/bin/bash

echo start

rm -r "VS_${1}"
#mkdir "VS_${1}"
#cd "VS_${1}"
pwd

samplename="/scratch/duhaimem_root/duhaimem/shared_data/VSTE/TestingSet/metagenomic_testing_set_${1}.fna"
threads=$2
path="/scratch/duhaimem_root/duhaimem/shared_data/VSTE/RunningToolsOnTestingSet/VS/VS_${1}"
vsdirectory="/nfs/turbo/lsa-duhaimem/software/virsorter-data"

/nfs/turbo/cee-kwigg/hegartyb/VirSorter/wrapper_phage_contigs_sorter_iPlant.pl \
        -f $samplename --db 1 \
        --wdir $path \
        --data-dir $vsdirectory \
        --ncpu $threads

echo done
