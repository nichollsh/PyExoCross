#!/bin/bash

# --------------------------
# CHANGE THESE

#SBATCH --mem=450GB
#SBATCH --time=02-22:58:00
#SBATCH --cpus-per-task=22
#SBATCH --job-name=xsec_SO2
#SBATCH --export=ALL
#SBATCH --output=slurm-%x-%j-stdout.log
#SBATCH --error=slurm-%x-%j-stderr.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=gpu
#SBATCH --gpus-per-node=a100:1

# CONFIG="input/NH3_ExoMol_CoYuTe_xsec.inp"
# CONFIG="input/CH4_ExoMol_MM_xsec.inp"
CONFIG="input/SO2_ExoMol_ExoAmes_xsec.inp"
RESUME="" # or "-r" to resume from checkpoint

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
echo "Resume:            $RESUME"

echo " "
srun python run.py -p $CONFIG $RESUME
echo " "

echo "Done"
