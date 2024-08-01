#!/bin/bash

#A typical run takes couple of hours but may be much longer
#SBATCH --job-name=array
#SBATCH --nodes=1
#SBATCH --time=02:00:00
#SBATCH --partition=gpu
#SBATCH --gpus=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --gpus-per-task=1

#log files:
#SBATCH -e logs/run_multimer_jobs_%A_%a_err.txt
#SBATCH -o logs/run_multimer_jobs_%A_%a_out.txt

TMP_DIR="/scratch-shared/tmp.KJVsRs6W2S"

TMP_OUTPUT="$TMP_DIR/output_ba2_100"

echo "Shared scratch space path: $TMP_DIR"

NODE_NAME=$(hostname)
echo "Node name: $NODE_NAME"

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
    --output_path=$TMP_OUTPUT \
    --data_dir=/projects/2/managed_datasets/AlphaFold \
    --protein_lists=baits.txt,candidates_ba2.txt \
    --monomer_objects_dir=./outputs/ba2_100 \
    --job_index=$SLURM_ARRAY_TASK_ID \
    --remove_result_pickles=False \
    

