#!/bin/bash

###############################
#                             #
#  1) Job Submission Options  #
#                             #
###############################

# Name
#SBATCH --job-name=startingruntoolstestingset

# Resources
# For MPI, increase ntasks-per-node
# For multithreading, increase cpus-per-task
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=01:00:00

# Account
#SBATCH --account=kwigg1
#SBATCH --partition=standard

# Logs
#SBATCH --mail-user=hegartyb@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# Environment
##SBATCH --export=ALL

#  Show list of CPUs you ran on
echo $SLURM_JOB_NODELIST

echo start

sbatch /scratch/duhaimem_root/duhaimem/shared_data/VSTE/RunningToolsOnTestingSet/Scripts/checkv_slurm.sh 
sbatch /scratch/duhaimem_root/duhaimem/shared_data/VSTE/RunningToolsOnTestingSet/Scripts/dvf_slurm.sh 
sbatch /scratch/duhaimem_root/duhaimem/shared_data/VSTE/RunningToolsOnTestingSet/Scripts/kaiju_slurm.sh 
sbatch /scratch/duhaimem_root/duhaimem/shared_data/VSTE/RunningToolsOnTestingSet/Scripts/vibrant_slurm.sh 
sbatch /scratch/duhaimem_root/duhaimem/shared_data/VSTE/RunningToolsOnTestingSet/Scripts/vs_slurm.sh 
sbatch /scratch/duhaimem_root/duhaimem/shared_data/VSTE/RunningToolsOnTestingSet/Scripts/vs2_slurm.sh

echo done
