#!/bin/bash

#A typical run takes couple of hours but may be much longer
#SBATCH --nodes=1
#SBATCH --time=03:00:00
#SBATCH --partition=genoa
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18

#SBATCH -e logs/run_feature_jobs_%A_%a_err.txt
#SBATCH -o logs/run_feature_jobs_%A_%a_out.txt

# module load 2022
# module load HH-suite/3.3.0-gompi-2022a
# module load HMMER/3.3.2-gompi-2022a
# module load Anaconda3/2022.05
module load 2023
module load Anaconda3/2023.07-2
source activate AlphaPulldown

create_individual_features.py \
  --fasta_paths=baits.fasta,candidates_ba2.fasta \
  --data_dir=/projects/2/managed_datasets/AlphaFold \
  --save_msa_files=False \
  --output_dir=./outputs/ba2_100 \
  --use_precomputed_msas=True \
  --max_template_date=2050-01-01 \
  --skip_existing=True \
  --seq_index=$SLURM_ARRAY_TASK_ID