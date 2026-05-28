#!/bin/bash

# LEAVE
#SBATCH --export=ALL
#SBATCH --job-name=pyexocross
#SBATCH --output=slurm-%x-%j-stdout.log
#SBATCH --error=slurm-%x-%j-stderr.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1

# CHANGE
#SBATCH --gpus-per-node=1
#SBATCH --mem=64GB
#SBATCH --time=00-03:00:00
#SBATCH --cpus-per-task=30
CONFIG="input/H2O_ExoMol_xsec.inp"


echo "Home dir:          $HOME"
source $HOME/.bashrc
conda activate exocross

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
