#!/bin/bash
#SBATCH --job-name=my_job
#SBATCH --partition rome
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=62
#SBATCH --time=15:00:00

#  

module load 2022
module load AlphaFold/2.3.1-foss-2022a-CUDA-11.7.0
module load parallel/20220722-GCCcore-11.3.0

fasta_input_dir="./inputs/ba2_100/"
output_base_dir="./outputs/ba2_100/"

# fasta_input_dir="outputs/msa/fastas/sequence_datasets/5/batch1"
# output_base_dir="outputs/msa/controls/batch1/"

export fasta_input_dir output_base_dir

run_commands() {
    fasta_file=$1
    fasta_filename=$(basename "$fasta_file")
    fasta_basename="${fasta_filename%.*}"

    # Create output directory for the fasta file
    output_dir="${output_base_dir}${fasta_basename}"
    mkdir -p "$output_dir"

    # Run hhblits
    hhblits -i "$fasta_file" -o "${output_dir}/pdb_hits.hhr" -cpu 8 -oa3m "${output_dir}/bfd_uniref_hits.a3m" -n 3 -e 0.001 -maxseq 1000000 -realign_max 100000 -maxfilt 100000 -min_prefilter_hits 1000 -d /projects/2/managed_datasets/AlphaFold/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt -d /projects/2/managed_datasets/AlphaFold/uniref30/UniRef30_2021_03

    # Run jackhmmer for uniref90
    jackhmmer -A "${output_dir}/uniref90_hits.sto" --noali --F1 0.0005 --F2 5e-05 --F3 5e-07 --incE 0.0001 -E 0.0001 --cpu 8 -N 1 "$fasta_file" /projects/2/managed_datasets/AlphaFold/uniref90/uniref90.fasta

    # Run jackhmmer for mgnify
    jackhmmer -A "${output_dir}/mgnify_output.sto" --noali --F1 0.0005 --F2 5e-05 --F3 5e-07 --incE 0.0001 -E 0.0001 --cpu 8 -N 1 "$fasta_file" /projects/2/managed_datasets/AlphaFold/mgnify/mgy_clusters_2022_05.fa
}

export -f run_commands

find "${fasta_input_dir}" -name "*.fasta" | parallel -j 12 run_commands {}