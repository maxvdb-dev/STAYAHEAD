#!/bin/bash

#A typical run takes couple of hours but may be much longer
#SBATCH --job-name=array
#SBATCH --nodes=1
#SBATCH --time=00:45:00
#SBATCH --partition=gpu
#SBATCH --gpus=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18

#log files:
#SBATCH -e logs/run_multimer_jobs_%A_%a_err.txt
#SBATCH -o logs/run_multimer_jobs_%A_%a_out.txt

module load 2023
module load Anaconda3/2023.07-2
module load 2022
module load CUDA/11.8.0
module load cuDNN/8.6.0.163-CUDA-11.8.0
source activate AlphaPulldown

MAXRAM=$(echo `ulimit -m` '/ 1024.0'|bc)
GPUMEM=`nvidia-smi --query-gpu=memory.total --format=csv,noheader,nounits|tail -1`
export XLA_PYTHON_CLIENT_MEM_FRACTION=`echo "scale=3;$MAXRAM / $GPUMEM"|bc`
export TF_FORCE_UNIFIED_MEMORY='1'

run_multimer_jobs.py --mode=pulldown \
    --num_cycle=3 \
    --num_predictions_per_model=1 \
    --output_path=./outputs/50_1-3 \
    --data_dir=/projects/2/managed_datasets/AlphaFold \
    --protein_lists=baits.txt,candidates_50_1-3.txt \
    --monomer_objects_dir=./outputs/1step/feat50_1 \
    --job_index=$SLURM_ARRAY_TASK_ID \
    --remove_result_pickles=True \