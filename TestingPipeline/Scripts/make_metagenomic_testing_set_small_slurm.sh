#!/bin/bash

###############################
#                             #
#  1) Job Submission Options  #
#                             #
###############################

# Name
#SBATCH --job-name=testingset

# Resources
# For MPI, increase ntasks-per-node
# For multithreading, increase cpus-per-task
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1000mb
#SBATCH --time=01-00:00:00

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

echo "started"

for i in {1..10}
do
    echo $i
    python make_metagenomic_testing_set_small.py $i
done

echo "done"
