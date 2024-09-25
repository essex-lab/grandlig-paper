#!/bin/bash
#
# Template script for SLURM jobs
# Note that "#SBATCH" is not a comment, but "# SBATCH" is (note the space after "#")
# For the online manual for sbatch, type "man sbatch" in the terminal

# Job name:
#----------------
#SBATCH --job-name="EQ"
#SBATCH --array=0-19
#----------------

# Account name (semi-optional):
#----------------
# SBATCH --account=
#----------------

# Partition name:
#----------------
#SBATCH --partition=gtx1080
#----------------
# Specifying resources needed for run:
#
#--------------
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --time=06:00:00
##SBATCH --mem=1gb
#--------------

# ---------------- ARRAY JOB SETTINGS ------------------ #
LAMBDA=$SLURM_ARRAY_TASK_ID
# -------------------------------------------------------#

#--------------
# Define some simpler local evironment variables:
nprocs=$SLURM_NTASKS
nnodes=$SLURM_JOB_NUM_NODES

#Function to call to run the actual code
slurm_startjob(){
#----------------- Actual calculation command goes here: ---------------------------

source ~/.bashrc
conda activate OMM8.0

# Set some environment variables
FREE_ENERGY=$SLURM_SUBMIT_DIR
echo "Free energy home directory set to $FREE_ENERGY"
echo $SLURM_NTASKS

# A new directory will be created for each value of lambda and
# at each step in the workflow for maximum organization.

mkdir lambda_$LAMBDA
cd lambda_$LAMBDA


SCRIPT_FOLDER=~/gcmc_scripts/FragmentScreening

/scratch/wp1g16/miniconda3/envs/OMM8.0/bin/python $SCRIPT_FOLDER/calcMuEx.py -p ../../cosolv_box.pdb -x ../../*.xml -lam $LAMBDA

echo "Ending. Job completed for lambda = $LAMBDA"


cd $FREE_ENERGY


echo Job Done
}

# Function to echo informational output
slurm_info_out(){
# Informational output
echo "=================================== SLURM JOB ==================================="
echo
echo "The job will be started on the following node(s):"
echo $SLURM_JOB_NODELIST
echo
echo "Slurm User:         $SLURM_JOB_USER"
echo "Run Directory:      $(pwd)"
echo "Job ID:             ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "Job Name:           $SLURM_JOB_NAME"
echo "Partition:          $SLURM_JOB_PARTITION"
echo "Number of nodes:    $SLURM_JOB_NUM_NODES"
echo "Number of tasks:    $SLURM_NTASKS"
echo "Submitted From:     $SLURM_SUBMIT_HOST:$SLURM_SUBMIT_DIR"
echo "=================================== SLURM JOB ==================================="
echo
echo "--- SLURM job-script output ---"
}

slurm_info_out

slurm_startjob
