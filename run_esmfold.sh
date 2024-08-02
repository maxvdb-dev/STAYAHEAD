#!/bin/bash
#SBATCH --job-name=esmfold_job
#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH --partition=gpu
#SBATCH --gpus=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH -e esmfold_%j_%t.out
#SBATCH -e esmfold_%j_%t.err

# Activate environment and load necessary modules
source ./venv/bin/activate
module load 2022
module load cuDNN/8.6.0.163-CUDA-11.8.0
cd ./esmfold/

# Define command arguments
cmd_args="--fastas_folder ./inputs/ds1/1273 \
--output_folder ./outputs/ds1/1273"

# Run ESMFold
python esmfold.py ${cmd_args}

# /usr/bin/time -v python esmfold.py ${cmd_args} 2>&1 | tee esmfold_time.log
