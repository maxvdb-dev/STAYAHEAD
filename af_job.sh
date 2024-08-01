#!/bin/bash
#SBATCH --job-name=my_job
#SBATCH --nodes=1
#SBATCH --time=15:00:00
#SBATCH --partition gpu
#SBATCH --gpus 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 18

module load 2022
module load AlphaFold/2.3.1-foss-2022a-CUDA-11.7.0

data_root=/projects/2/managed_datasets/AlphaFold
project_root=./

fasta_files=("${project_root}outputs/msa/fastas/missings_seqs_fastas/"*.fasta)

cmd_args="--fasta_paths {}
    --max_template_date 2023-03-20
    --data_dir ${data_root}
    --db_preset full_dbs
    --output_dir ${project_root}outputs/msa/missing_seqs_msas/{/.}/{/.}
    --use_precomputed_msas" ::: "${fasta_files[@]}"

alphafold ${cmd_args}
