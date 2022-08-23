#!/bin/bash

###############################
#                             #
#  1) Job Submission Options  #
#                             #
###############################

# Name
#SBATCH --job-name=kaiju

# Resources
# For MPI, increase ntasks-per-node
# For multithreading, increase cpus-per-task
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=20G
#SBATCH --time=5-00:00:00

# Account
#SBATCH --account=kwigg1
#SBATCH --partition=standard

# Logs
#SBATCH --mail-user=hegartyb@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# Environment
##SBATCH --export=ALL

#SBATCH --array=1-10

#  Show list of CPUs you ran on
echo $SLURM_JOB_NODELIST

echo start

echo "activate conda environment"
source activate  /home/hegartyb/miniconda3/envs/cluster_kaiju

echo "run dvf on all testing sets"
mkdir ../Kaiju
cd ../Kaiju
bash ../Scripts/kaiju.sh ${SLURM_ARRAY_TASK_ID} ${SLURM_CPUS_ON_NODE}

echo done
