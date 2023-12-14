





##############################
# 1 - Load libraries
##############################

library(ggplot2)
library(stringr)
library(plyr)
library(reshape2)
library(viridis)
library(tidyr)
library(dplyr)
library(readr)
library(data.table)

##############################
# 2 - Set parameters
##############################

option_list = list(
  make_option(c("-i", "--infile"), type="character", default="./merged_viral_feature_outputs.csv", 
              help="input path/filename from 01-build-feature-table.R [default = %default]", metavar="character"),
  
  make_option(c("-a", "--virsorter2"), type="character", default=FALSE, 
              help="specify if <tool_name> was included in <infile> generation [default = %default]", metavar="character"),
  
  make_option(c("-b", "--virsorter"), type="character", default=FALSE, 
              help="specify if <tool_name> was included in <infile> generation [default = %default]", metavar="character"),
  
  make_option(c("-c", "--vibrant"), type="character", default=FALSE, 
              help="specify if <tool_name> was included in <infile> generation [default = %default]", metavar="character"),
  
  make_option(c("-d", "--deepvirfinder"), type="character", default=FALSE, 
              help="specify if <tool_name> was included in <infile> generation [default = %default]", metavar="character"),
  
  make_option(c("-t", "--tuning_viral"), type="character", default=FALSE, 
              help="specify if <tuning_viral> is wished to be used. Requires virsorter2, checkv, and kaiju. [default = %default]", metavar="character"),
  
  make_option(c("-n", "--tuning_not_viral"), type="character", default=FALSE, 
              help="specify if <tuning_viral> is wished to be used. Requires CheckV. [default = %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);



##############################
# 2 - compute viral score
##############################

# This section defines a viralness score "keep_score" based on the tool classifications. 
# A final keep_score above 1 indicates we will keep that sequence and call it viral.
# 
# Rule
# subrule
# ...
# subrule
# 
# VIBRANT
# Quality == "Complete Circular": +1
# Quality == "High Quality Draft": +1
# Quality == "Medium Quality Draft": +1
# Quality == "Low Quality Draft" & provirus == TRUE: +0.5
# 
# Virsorter
# category ==  1 or 4: +1
# category == 2, 3, 5, or 6: +0.5
# 
# Virsorter2
# max_score >= 0.5: +0.5
# max_score >= 0.95: +0.5 (additive to the previous subrule)
# 
# DeepVirFinder:
#   Score >= 0.7 & p-value <= 0.05 & checkv_length < 20kb: +0.5
# Score >= 0.9 & p-value <= 0.05 & checkv_length < 20kb: +0.5 (additive to the previous subrule)
# 
# Tuning Addition:
#   Kaiju_viral == "Viruses" & kaiju_match_ratio >= 0.3: +1
# If %unknown >= 75 and checkv_length < 50000: +0.5
# If checkv_%viral >= 50 or vs2_%viral >= 50%: +0.5
# Hallmark >= 3: +1
# 
# Tuning - Not Viral:
#   viral_genes == 0 & host_genes >= 1: -1
# If host_genes >50 and checkv_provirus == False: -1 
# If 3*viral_genes <= host_genes and  checkv_provirus == False: -1
# If length > 50,000 and hallmark <=1: -1
# 
# 
# This script produces visualizations of these combined viral scorings and
# includes ecological metrics like alpha diversity.
# 
# Users can decide which combination is appropriate for them and only need to run
# the chunks for the appropriate tools.

compute_viral_score <- function(input_seqs,
                                include_vibrant=F, 
                                include_virsorter=F,
                                include_virsorter2=F,
                                include_deepvirfinder=F,
                                include_tuning_viral=F,
                                tv_1=T, tv_2=T, tv_3=T, tv_4=T,
                                include_tuning_not_viral=F,
                                ntv_1=T, ntv_2=T, ntv_3=T, ntv_4=T) {
  
  keep_score <- rep(0, nrow(input_seqs))
  if (include_vibrant) {
    keep_score[input_seqs$vibrant_quality=="high quality draft"] <- keep_score[input_seqs$vibrant_quality=="high quality draft"] + 1
    keep_score[input_seqs$vibrant_quality=="medium quality draft"] <- keep_score[input_seqs$vibrant_quality=="medium quality draft"] + 1
    keep_score[input_seqs$vibrant_quality=="low quality draft" & input_seqs$vibrant_prophage] <- keep_score[input_seqs$vibrant_quality=="low quality draft" & input_seqs$vibrant_prophage] + 0.5
  }
  
  if (include_virsorter2) {
    keep_score[input_seqs$max_score>=0.50] <- keep_score[input_seqs$max_score>=0.50] + 0.5
    keep_score[input_seqs$max_score>=0.95] <- keep_score[input_seqs$max_score>=0.95] + 0.5  
  }
  
  if (include_virsorter) {
    keep_score[input_seqs$category==1] <- keep_score[input_seqs$category==1] + 1
    keep_score[input_seqs$category==2] <- keep_score[input_seqs$category==2] + 0.5
    keep_score[input_seqs$category==3] <- keep_score[input_seqs$category==3] + 0.5
    keep_score[input_seqs$category==4] <- keep_score[input_seqs$category==4] + 1
    keep_score[input_seqs$category==5] <- keep_score[input_seqs$category==5] + 0.5 
    keep_score[input_seqs$category==6] <- keep_score[input_seqs$category==6] + 0.5 
  }
  
  if (include_deepvirfinder) {
    keep_score[(input_seqs$score>=0.9 & input_seqs$pvalue<=0.05) & input_seqs$checkv_length<20000] <- keep_score[(input_seqs$score>=0.9 & input_seqs$pvalue<=0.05) & input_seqs$checkv_length<20000] + 0.5
    keep_score[(input_seqs$score>=0.7 & input_seqs$pvalue<=0.05) & input_seqs$checkv_length<20000] <- keep_score[(input_seqs$score>=0.7 & input_seqs$pvalue<=0.05) & input_seqs$checkv_length<20000] + 0.5
  }
  
  if (include_tuning_viral) {
    keep_score <- tuning_viral(current_score=keep_score,
                               ins=input_seqs,
                               include_tv_1=tv_1,
                               include_tv_2=tv_2,
                               include_tv_3=tv_3,
                               include_tv_4=tv_4)
  }
  
  if (include_tuning_not_viral) {
    # tuning removal 
    keep_score <- tuning_not_viral(current_score=keep_score,
                                   ins=input_seqs,
                                   include_tnv_1=ntv_1,
                                   include_tnv_2=ntv_2,
                                   include_tnv_3=ntv_3,
                                   include_tnv_4=ntv_4)
    
  }
  
  return(keep_score)
  
}
tuning_viral <- function(current_score, ins, 
                         include_tv_1=T,
                         include_tv_2=T,
                         include_tv_3=T,
                         include_tv_4=T) {
  # tuning addition
  if (include_tv_1) {
    current_score[ins$hallmark>=3] <- current_score[ins$hallmark>=3] + 1
  }
  if (include_tv_2) {
    current_score[ins$Kaiju_Viral=="Viruses" & ins$kaiju_match_ratio>=0.3] <- current_score[ins$Kaiju_Viral=="Viruses" & ins$kaiju_match_ratio>=0.3] + 1
  }
  if (include_tv_3) {
    current_score[ins$percent_unknown>=75 & ins$checkv_length<50000] <- current_score[ins$percent_unknown>=75 & ins$checkv_length<50000] + 0.5
  }
  if (include_tv_4) {
    current_score[ins$viral>=50 | ins$percent_viral>=50] <- current_score[ins$viral>=50 | ins$percent_viral>=50] + 0.5 
  }
  return(current_score)
}
tuning_not_viral <- function(current_score, ins, 
                             include_tnv_1=T,
                             include_tnv_2=T,
                             include_tnv_3=T,
                             include_tnv_4=T) {
  # tuning removal
  if (include_tnv_1) {
    current_score[(ins$checkv_host_genes>50) & !ins$checkv_provirus] <- -3
  }
  
  if (include_tnv_2) {
    current_score[ins$checkv_viral_genes==0 & ins$checkv_host_genes>=1] <- current_score[ins$checkv_viral_genes==0 & ins$checkv_host_genes>=1] - 1
  }
  if (include_tnv_3) {
    current_score[((ins$checkv_viral_genes*3) <= ins$checkv_host_genes) & !ins$checkv_provirus] <- current_score[((ins$checkv_viral_genes*3) <= ins$checkv_host_genes) & !ins$checkv_provirus] - 1
  }
  
  if (include_tnv_4) {
    current_score[ins$checkv_length>500000 & ins$hallmark<=1] <- -3
  }
  
  return(current_score)
}

###################################
# 2 - apply viral score function
###################################

# The next code chunk is where you can create your own combination of tools,
# or use one of the presets provided. Simply change each include_* rule to True
# to include the rule, and False to exclude it.


command_line_params <- paste(names(opt)[names(opt) != "infile" & sapply(opt, function(x) x != FALSE)], collapse="__")

viruses[[command_line_params]] <- compute_viral_score(viruses,
                                                      include_vibrant=opt$vibrant, 
                                                      include_virsorter=opt$virsorter,
                                                      include_virsorter2=opt$virsorter2,
                                                      include_deepvirfinder=opt$deepvirfinder,
                                                      include_tuning_viral=opt$tuning_viral,
                                                      include_tuning_not_viral = opt$tuning_not_viral)
viruses$keep_score_all <- compute_viral_score(viruses, include_deepvirfinder = T,
                                              include_vibrant = T,
                                              include_virsorter2 = T,
                                              include_tuning_viral = T,
                                              include_tuning_not_viral = T,
                                              include_virsorter = T)
viruses$keep_score_high_recall <- compute_viral_score(viruses, include_deepvirfinder = T,
                                                      include_vibrant = T,
                                                      include_virsorter2 = T,
                                                      include_tuning_viral = T,
                                                      include_tuning_not_viral = F,
                                                      include_virsorter = T)
viruses$keep_score_high_MCC <- compute_viral_score(viruses, include_deepvirfinder = F,
                                                   include_vibrant = F,
                                                   include_virsorter2 = T,
                                                   include_tuning_viral = F,
                                                   include_tuning_not_viral = T,
                                                   include_virsorter = F)
viruses$keep_score_high_precision <- compute_viral_score(viruses, include_deepvirfinder = F,
                                                         include_vibrant = F,
                                                         include_virsorter2 = T,
                                                         include_tuning_viral = F,
                                                         include_tuning_not_viral = T,
                                                         include_virsorter = F)

