#!/bin/bash

###############################
#                             #
#  1) Job Submission Options  #
#                             #
###############################

# Name
#SBATCH --job-name=overarching

# Resources
# For MPI, increase ntasks-per-node
# For multithreading, increase cpus-per-task
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=300mb
#SBATCH --time=04-00:00:00

# Account
#SBATCH --account=kwigg1
#SBATCH --partition=standard

# Logs
#SBATCH --mail-user=hegartyb@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=Logs/%x-%j.out

# Environment
##SBATCH --export=ALL

source /etc/profile.d/http_proxy.sh

# List compute nodes allocated to the job
if [[ $SLURM_JOB_NODELIST ]] ; then
    echo "Running on"
    scontrol show hostnames $SLURM_JOB_NODELIST
    echo -e "\n"
fi



#####################
#                   #
#  2) Job Commands  #
#                   #
#####################

# Initiating snakemake and running workflow in cluster mode
snakemake --profile /scratch/kwigg_root/kwigg/hegartyb/SnakemakeAssemblies3000/CompetitiveMapping/Config --latency-wait 60 --use-conda --conda-prefix /home/hegartyb/miniconda3/envs/ --snakefile Snakefile-overarching --keep-going --rerun-incomplete




