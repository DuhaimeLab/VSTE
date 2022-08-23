#!/bin/bash

echo start

rm -r "Kaiju_${1}"
mkdir "Kaiju_${1}"
cd "Kaiju_${1}"
pwd

samplename="/scratch/duhaimem_root/duhaimem/shared_data/VSTE/TestingSet/metagenomic_testing_set_${1}.fna"
threads=$2
input_ndmp="/scratch/duhaimem_root/duhaimem0/shared_data/kaijudb/nr_euk/nodes.dmp"
input_names="/scratch/duhaimem_root/duhaimem0/shared_data/kaijudb/nr_euk/names.dmp"
input_vmi="/nfs/turbo/cee-kwigg/hegartyb/Kaiju/Kaijudb-refseq/kaiju_db_refseq.fmi"
outname="${1}.nreuk.kaiju.out"
outname_taxa="${1}.nreuk.kaiju.names.out"

/nfs/turbo/lsa-duhaimem/Lake_Michigan_JGI_viral_metaGs/kaiju/bin/kaiju -z $2 -t $input_ndmp -f $input_vmi -i $samplename -o $outname -v
/nfs/turbo/lsa-duhaimem/Lake_Michigan_JGI_viral_metaGs/kaiju/bin/kaiju-addTaxonNames -t $input_ndmp -n $input_names -i $outname -o $outname_taxa -p

echo done
