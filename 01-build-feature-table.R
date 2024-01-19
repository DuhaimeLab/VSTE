#!/usr/bin/env Rscript

#############################################################
## Takes outputs from viral identification tools and merges
## them into one file, "merged_viral_id_outputs.csv", which is the input
## for identify_viruses.R
##
## NOTE: This script is a blueprint for cleaning up the feature table. Your table may
## be different from ours and will require modifying the script for the data cleaning process
## I have tried my bes to make it as user friendly as possible, but will require some modifications
## depending on how your files are named and which parameters and versions 
## were used for each viral ID output, since these can influence the filenames.
## These changes to the data wrangling process can be made in the 2 - load data section of this script.
##
## Please reach out to James Riddell (riddell.26@buckeyemail.osu.edu) or
## Bridget Hegarty (beh53@case.edu) regarding any issues, or open an issue on github.
#############################################################


#############################################################
# Please read the below documentation to ensure your input
# files are in the correct format before merging.
#############################################################

# 1) For this script's inputs, the user must merge each viral identification
# tool output with the following specified filenames into one tab-separated file.
# Replace the filenames with their paths if needed into the variables in the
# 'inputs' chunk.
# 
#     General output filenames used from each tool:
#         Checkv: quality_summary.tsv
#         VIBRANT: VIBRANT_genome_quality_${assembly}.tsv
#         Virsorter: VIRSorter_global-phage-signal.csv
#         Virsorter2: final-viral-score.tsv
#         DeepVirFinder: ${assembly}.fasta_gt2500bp_dvfpred.txt
#         Kaiju: ${assembly}.kaiju.names.out
# 
# 2) Make sure all fasta headers (assembly, contig) are consistent within and
# across tools. Chunks for each tool do contain some lines for cleaning up these features, but due to their variability it will be the user's responsibility to
# make sure they match across tools.
# 
# 3) Check each chunk and ensure all columns are accounted for. For example,
# if you don't have an "assembly" column and called it something else, please either
# remove it from this script or modify your input files.
# 
# 4) This script is designed for contigs > 3000 basepairs. It can be modified
# to be higher or lower, but going lower will greatly increase the size of the
# dataframe and memory usage.


##############################
# 1 - Load libraries
##############################
library(ggplot2)
library(plyr)
library(reshape2)
library(viridis)
library(tidyr)
library(dplyr)
library(readr)
library(data.table)
library(optparse)


##############################
# 1 - Set parameters
##############################

option_list = list(
  make_option(c("-o", "--outfile"), type="character", default="./merged_viral_feature_outputs.csv", 
              help="output filename/path [default = %default]", metavar="character"),
  
  make_option(c("-l", "--length"), type="integer", default=3000, 
              help="length cutoff in basepairs [default= %default]", metavar="integer"),
  
  make_option(c("-a", "--virsorter2"), type="character", default="skip", 
              help="path to final-viral-score.tsv (including filename) [default = %default", metavar="character"),
  
  make_option(c("-b", "--virsorter"), type="character", default="skip", 
              help="path to VIRSorter_global-phage-signal.csv (including filename) [default = %default", metavar="character"),
  
  make_option(c("-c", "--vibrant"), type="character", default="skip", 
              help="path to VIBRANT_genome_quality_${samplename}.tsv (including filename) [default = %default", metavar="character"),
  
  make_option(c("-d", "--deepvirfinder"), type="character", default="skip", 
              help="path to ${samplename}.fasta_gt${cutoff}bp_dvfpred.txt (including filename) [default = %default", metavar="character"),
  
  make_option(c("-e", "--checkv"), type="character", default="skip", 
              help="path to quality_summary.tsv (including filename) [default = %default", metavar="character"),
  
  make_option(c("-f", "--kaiju"), type="character", default="skip", 
              help="path to ${samplename}.kaiju.names.out (including filename) [default = %default", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Virsorter2
vs2_path <- opt$virsorter2

# Virsorter
vs_path <- opt$virsorter

# VIBRANT
vb_path <- opt$vibrant

# DeepVirFinder
dvf_path <- opt$deepvirfinder

# cv_c
cv_path <- opt$checkv

# Kaiju
kj_path <- opt$kaiju

# length cutoff (unit: basepairs)
BP_CUTOFF <- opt$length

##############################
# 2 - load data
##############################

# add dataframes if they are detected
df_counter <- 0
df_c <- list()

if (cv_path != "skip") {
  cv_c <- fread(cv_path, 
                  sep="\t",
                  header = T, 
                  select = c(
                    'contig_id',
                    'assembly',
                    'provirus',
                    'completeness',
                    'contamination',
                    'viral_genes',
                    'host_genes',
                    'gene_count',
                    'contig_length',
                    'checkv_quality'
                  )
  ) %>% 
    rename(
      contig = contig_id,
      checkv_provirus = provirus,
      checkv_completeness = completeness,
      checkv_contamination = contamination,
      checkv_viral_genes = viral_genes,
      checkv_host_genes = host_genes,
      checkv_total_genes = gene_count,
      checkv_length = contig_length
    )
  cv_c$checkv_provirus <- cv_c$checkv_provirus == 'Yes'
  cv_c$method = 'checkv'
  cv_c$contig <- sub("\\.", "_", cv_c$contig)
  cv_c$uniq_contig <- paste(cv_c$assembly, cv_c$contig, sep = "--")
  cv_c <- cv_c[!duplicated(cv_c$uniq_contig),]
  
  df_c[[df_counter+1]] <- cv_c
  df_counter <- df_counter + 1
} else {
  print("skipping CheckV, no CheckV output was detected.")
}




if (vb_path != "skip") {
  vb_c <- fread(vb_path,
                header = T,
                sep = "\t",
                select = c(
                  'scaffold',
                  'assembly',
                  'type',
                  'Quality'
                )
  ) %>%
    rename(
      contig = scaffold,
      vibrant_quality = Quality
    )
  vb_c$method <- "vibrant"
  vb_c$vibrant_prophage <- FALSE
  vb_c$vibrant_prophage[grep("_fragment_", vb_c$contig)] <- TRUE
  vb_c$contig <- gsub("_fragment_.*", "", vb_c$contig)
  vb_c$contig <- sub("\\.", "_", vb_c$contig)
  vb_c$uniq_contig <- paste(vb_c$assembly, vb_c$contig, sep = "--")
  vb_c <- vb_c[!duplicated(vb_c$uniq_contig),]
  
  df_c[[df_counter+1]] <- vb_c
  df_counter <- df_counter + 1
} else {
  vb_c <- data.frame()
  print("skipping VIBRANT, no VIBRANT output was detected.")
}


if (dvf_path != "skip") {
  dvf_c <- fread(dvf_path,
                 header = T,
                 sep = "\t",
                 select = c(
                   'assembly',
                   'name',
                   'score',
                   'pvalue'
                 )
  ) %>% 
    rename(
      contig = name,
    )
  dvf_c$contig <- sub("\\.", "_", dvf_c$contig)
  dvf_c$bh_pvalue <- p.adjust(dvf_c$pvalue, method = "BH")
  dvf_c$uniq_contig <- paste(dvf_c$assembly, dvf_c$contig, sep = "--")
  dvf_c <- dvf_c[!duplicated(dvf_c$uniq_contig),]
  
  df_c[[df_counter+1]] <- dvf_c
  df_counter <- df_counter + 1
} else {
  dvf_c <- data.frame()
  print("skipping DeepVirFinder, no DeepVirFinder output was detected.")
}


if (vs_path != "skip") {
  vs_c <- fread(vs_path,
                select = c(
                  'assembly',
                  '## Contig_id',
                  'Category'
                ),
  ) %>% rename(uniq_contig = '## Contig_id') %>% drop_na(uniq_contig)
  vs_c$uniq_contig <- sub("VIRSorter_", "", vs_c$uniq_contig)
  vs_c$uniq_contig <- sub("-circular", "", vs_c$uniq_contig)
  vs_c$uniq_contig <- sub("\\.", "_", vs_c$uniq_contig)
  
  vs_c <- vs_c[!duplicated(vs_c$uniq_contig),]
  vs_c <- vs_c %>% drop_na(uniq_contig) 
  
  df_c[[df_counter+1]] <- vs_c
  df_counter <- df_counter + 1
} else {
  vs_c <- data.frame()
  print("skipping VirSorter, no VirSorter output was detected.")
}
 

if (vs2_path != "skip") {
  vs2_c <- fread(vs2_path,
                 header = T,
                 sep = '\t',
                 select = c(
                   'seqname',
                   'assembly',
                   'dsDNAphage',
                   'ssDNA',
                   'max_score',
                   'max_score_group',
                   'hallmark',
                   'viral',
                   'cellular'
                 )
  ) %>% 
    separate(
      col = seqname,
      into = c("contig", "type"), 
      sep = "\\|\\|",
      remove = T
    )
  vs2_c$contig <- sub("\\.", "_", vs2_c$contig)
  vs2_c$uniq_contig <- paste(vs2_c$assembly, vs2_c$contig, sep="--")
  vs2_c <- vs2_c[!duplicated(vs2_c$uniq_contig),]
  
  df_c[[df_counter+1]] <- vs2_c
  df_counter <- df_counter + 1
} else {
  vs2_c <- data.frame()
  print("skipping VirSorter2, no VirSorter2 output was detected.")
}


if (kj_path != "skip") {
  kj_c <- fread(kj_path,
                 header = T,
                 sep = '\t',
                 select = c(
                   'assembly_contig',
                   'length_or_score_of_best_match',
                   'taxonomy',
                   'assembly'
                 )
  ) %>% rename(contig = assembly_contig)
  
  kj_c <- separate(kj_c, col = taxonomy, into = c("Kaiju_Viral","Kingdom"), sep=";")
  kj_c$uniq_contig <- paste(kj_c$assembly, kj_c$contig, sep="--")
  kj_c <- kj_c[!duplicated(kj_c$uniq_contig),]
  kj_c$uniq_contig <- sub("\\.", "_", kj_c$uniq_contig)
  
  df_c[[df_counter+1]] <- kj_c
  df_counter <- df_counter + 1
} else {
  kj_c <- data.frame()
  print("skipping Kaiju, no Kaiju output was detected.")
}

##############################
# 3 - merge tool outputs
##############################

library(purrr)
viruses <- df_c %>% reduce(full_join, by=c("uniq_contig", "assembly"))

##############################
# 4 - filter by length cutoff
##############################

if (nrow(cv_c) != 0) {
viruses <- viruses %>% filter(checkv_length > BP_CUTOFF)
}

##############################
# 5 - check missing
##############################

# Use the next code chunk to check if any contigs are missing from checkv.
# If there are any rows present, these are located in other viral ID outputs,
# but not CheckV. Please go back and change the contig names so that they match
# the pattern in CheckV, or check your CheckV output to make sure it contains
# as many rows as headers in your fasta file (with the bp cutoff).

if (nrow(cv_c) != 0) {
  v_missing <- viruses[is.na(viruses$checkv_quality),]
} else {
  print("Cannot check for missing sequences in the merged viruses dataframe because no CheckV output was detected.")
}


###############################################################################
# 6 - Calculate percent viral/host/unknown/kaiju_match_ratio columns
#     These will be used in the classification step.
###############################################################################

if (nrow(cv_c) != 0) {
  viruses$percent_host <- viruses$checkv_host_genes/viruses$checkv_total_genes*100
  viruses$percent_viral <- viruses$checkv_viral_genes/viruses$checkv_total_genes*100
  viruses$percent_unknown <- 100-(viruses$checkv_host_genes+viruses$checkv_viral_genes)/viruses$checkv_total_genes*100
} else {
  print("Cannot calculate percent host, viral, and unknown features because no checkV output was detected.")
}

if (nrow(kj_c) != 0) {
  viruses$kaiju_match_ratio <- viruses$length_or_score_of_best_match/viruses$checkv_length
} else {
  print("Cannot calculate kaiju match ratio, no kaiju output was detected.")
}
###############################################################################
# 7 - remove NAs
# Make sure this step is not removing NAs that shouldn't be NA! 
# Check section 2 dataframes to make sure no columns are generating NAs
###############################################################################

# checkV
if (nrow(cv_c) != 0) {
  viruses$percent_viral[is.na(viruses$percent_viral)] <- 0
  viruses$percent_unknown[is.na(viruses$percent_unknown)] <- 0
  viruses$checkv_completeness[is.na(viruses$checkv_completeness)] <- 0
}

# VIBRANT
if (nrow(vb_c) != 0) {
  viruses$vibrant_quality[is.na(viruses$vibrant_quality)] <- 0
}

# DeepVirFinder
if (nrow(dvf_c) != 0) {
  viruses$score[is.na(viruses$score)] <- 0
  viruses$bh_pvalue[is.na(viruses$bh_pvalue)] <- 0
}

# VirSorter2
if (nrow(vs2_c) != 0) {
  viruses$viral[is.na(viruses$viral)] <- 0
  viruses$hallmark[is.na(viruses$hallmark)] <- 0
  viruses$max_score[is.na(viruses$max_score)] <- 0
}

#Virsorter
if (nrow(vs_c) != 0) {
  viruses$Category[is.na(viruses$Category)] <- 0
}

# Kaiju
if (nrow(kj_c) != 0) {
  viruses$Kaiju_Viral[is.na(viruses$Kaiju_Viral)] <- "unknown"
  viruses$kaiju_match_ratio[is.na(viruses$kaiju_match_ratio)] <- 0
}

##############################
# 5 - export
##############################

viruses_filename <- opt$outfile
write.csv(viruses, viruses_filename)
print(paste('Wrote merged viruses dataframe to ',viruses_filename, sep=""))

