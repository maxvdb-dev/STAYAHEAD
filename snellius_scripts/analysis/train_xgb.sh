#!/bin/bash
#SBATCH --job-name=xgboost
#SBATCH --nodes=1
#SBATCH --time=120:00:00
#SBATCH --partition gpu
#SBATCH --gpus 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 18

#SBATCH -e logs/train_xgb_job_1H_20B_err.txt
#SBATCH -o logs/train_xgb_job_1H_20B_out.txt

module load 2023
module load Python/3.11.3-GCCcore-12.3.0
module load scikit-learn/1.3.1-gfbf-2023a
module load matplotlib/3.7.2-gfbf-2023a

python full_analysis.py -i ./inputs/1H_20B -o ./outputs/1H_20B