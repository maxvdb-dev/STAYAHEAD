#!/bin/bash
#SBATCH --job-name=my_job
#SBATCH --partition rome
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=02:00:00

cd ~/project/alphafold

module load 2022
module load AlphaFold/2.3.1-foss-2022a-CUDA-11.7.0

outputs/msa/fastas/sequence_datasets/4/ace2

hhblits -i ./outputs/msa/fastas/6M0J.fasta -cpu 8 -oa3m ./outputs/msa/bfd_uniref_hits.a3m -n 3 -e 0.001 -maxseq 1000000 -realign_max 100000 -maxfilt 100000 -min_prefilter_hits 1000 -d /projects/2/managed_datasets/AlphaFold/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt -d /projects/2/managed_datasets/AlphaFold/uniref30/UniRef30_2021_03

jackhmmer -A ./outputs/msa/uniref90_hits.sto --noali --F1 0.0005 --F2 5e-05 --F3 5e-07 --incE 0.0001 -E 0.0001 --cpu 8 -N 1 ./outputs/msa/fastas/6M0J.fasta /projects/2/managed_datasets/AlphaFold/uniref90/uniref90.fasta

jackhmmer -A ./outputs/msa/mgnify_output.sto --noali --F1 0.0005 --F2 5e-05 --F3 5e-07 --incE 0.0001 -E 0.0001 --cpu 8 -N 1 ./outputs/msa/fastas/6M0J.fasta /projects/2/managed_datasets/AlphaFold/mgnify/mgy_clusters_2022_05.fa