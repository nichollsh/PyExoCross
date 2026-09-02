#!/bin/bash

# --------------------------
# CHANGE THESE

#SBATCH --mem=400GB
#SBATCH --time=02-12:00:00
#SBATCH --cpus-per-task=50
#SBATCH --job-name=xsec_OCS
#SBATCH --export=ALL
#SBATCH --output=slurm-%x-%j-stdout.log
#SBATCH --error=slurm-%x-%j-stderr.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1

CONFIG="input/OCS_ExoMol_OYT8_xsec.inp"

# ---------------------------
# Dispatch the job

echo "Home dir:          $HOME"
source $HOME/.bashrc
conda activate xsec
module load netCDF netCDF-Fortran

echo "Python:            $(which python)"
echo "TMPDIR:            $TMPDIR"
echo "Allocated CPUs:    $SLURM_CPUS_PER_TASK"
echo "Allocated GPUs:    $SLURM_GPUS_PER_NODE"
echo "Allocated node:    $SLURM_JOB_NODELIST"

CONFIG=$(realpath $CONFIG)
echo "Config file:       $CONFIG"

echo "Started at:        $(date)"
echo "Expected end:      $(date -d @$SLURM_JOB_END_TIME)    [EPOCH=$SLURM_JOB_END_TIME]"

echo " "
srun python run.py -p $CONFIG
echo " "

echo "Done"
