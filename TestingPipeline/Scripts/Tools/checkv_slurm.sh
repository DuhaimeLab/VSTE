#!/bin/bash

###############################
#                             #
#  1) Job Submission Options  #
#                             #
###############################

# Name
#SBATCH --job-name=checkv

# Resources
# For MPI, increase ntasks-per-node
# For multithreading, increase cpus-per-task
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=80G
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
source activate  /home/hegartyb/miniconda3/envs/checkV

echo "run checkV on all testing sets"
bash checkv.sh ${SLURM_ARRAY_TASK_ID} {SLURM_CPUS_ON_NODE}

echo done