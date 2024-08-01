#!/bin/bash
#SBATCH --job-name=my_job
#SBATCH --nodes=1
#SBATCH --time=04:00:00
#SBATCH --partition gpu
#SBATCH --gpus 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 18

module load 2022
module load AlphaFold/2.3.1-foss-2022a-CUDA-11.7.0


data_root=/projects/2/managed_datasets/AlphaFold
project_root=./

cmd_args="--fasta_paths ${project_root}inputs/SARS-CoV-ACE2.fasta
    --max_template_date 2023-03-20
    --db_preset full_dbs
    --data_dir ${data_root}
    --output_dir ${project_root}outputs"

alphafold ${cmd_args}