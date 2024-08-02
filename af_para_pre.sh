#!/bin/bash
#SBATCH --job-name=af_known_variants
#SBATCH --nodes=2
#SBATCH --time=06:00:00 
#SBATCH --partition=gpu
#SBATCH --gpus-per-node=4 
#SBATCH --tasks-per-node=4 
#SBATCH --cpus-per-task=18 

module load 2022
module load AlphaFold/2.3.1-foss-2022a-CUDA-11.7.0
module load parallel/20220722-GCCcore-11.3.0

data_root=/projects/2/managed_datasets/AlphaFold
project_root=./

# fasta_files=("${project_root}inputs/sequence_datasets/4/ace2/"*.fasta)
fasta_files=("${project_root}outputs/msa/fastas/1step/batch6/"*.fasta)
# msa_dirs=("${project_root}outputs/variants/"*)

parallel -j 8 alphafold --fasta_paths {} \
    --output_dir ${project_root}outputs/msa/1step/batch6/{/.}/{/.}/ \
    --db_preset full_dbs --data_dir ${data_root} \
    --max_template_date 2023-03-20 ::: "${fasta_files[@]}" \
    --use_precomputed_msas

# Use GNU parallel to run alphafold on each fasta file
# Replace 'N' with the number of jobs you want to run in parallel. It should not exceed the number of GPUs.
# parallel -j 26 alphafold --fasta_paths {} \
#     --output_dir ${project_root}outputs/ds4/{/.} \
#     --db_preset full_dbs --data_dir ${data_root} \
#     --max_template_date 2023-03-20 ::: "${fasta_files[@]}" \
#     --use_precomputed_msas ${project_root}outputs/variants/{/.}/ ::: "${fasta_files[@]}"

# Use GNU parallel to run alphafold on each fasta file
# Replace 'N' with the number of jobs you want to run in parallel. It should not exceed the number of GPUs.
# parallel -j 26 alphafold --use_precomputed_msas {} \
#     --output_dir ${project_root}outputs/ds4/{/} \
#     --db_preset full_dbs --data_dir ${data_root} \
#     --max_template_date 2023-03-20 ::: "${msa_dirs[@]}"