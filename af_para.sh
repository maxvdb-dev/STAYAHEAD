#!/bin/bash
#SBATCH --job-name=af_known_variants
#SBATCH --nodes=1
#SBATCH --time=08:00:00 
#SBATCH --partition=gpu
#SBATCH --gpus-per-node=4 
#SBATCH --tasks-per-node=4 
#SBATCH --cpus-per-task=18 

module load 2022
module load AlphaFold/2.3.1-foss-2022a-CUDA-11.7.0
module load parallel/20220722-GCCcore-11.3.0

data_root=/projects/2/managed_datasets/AlphaFold
project_root=./

# fasta_files=(
#     "${project_root}inputs/sequence_datasets/4/ace2/eta_variant.fasta"
#     "${project_root}inputs/sequence_datasets/4/ace2/gamma_variant.fasta"
#     "${project_root}inputs/sequence_datasets/4/ace2/iota_variant.fasta"
#     "${project_root}inputs/sequence_datasets/4/ace2/kappa_variant.fasta"
# )

fasta_files=("${project_root}inputs/variants"*.fasta)

# Use GNU parallel to run alphafold on each fasta file
# Replace 'N' with the number of jobs you want to run in parallel. It should not exceed the number of GPUs.
parallel -j 4 alphafold --fasta_paths {} \
    --output_dir ${project_root}outputs/variants/{/.} \
    --db_preset full_dbs --data_dir ${data_root} \
    --max_template_date 2023-03-20 ::: "${fasta_files[@]}"