#Unused Lines:

original_fasta="${project_root}inputs/10/sequences_10.fasta"

# Directory to store the split fasta files
split_fasta_dir="${project_root}inputs/split_fastas"
mkdir -p ${split_fasta_dir}

# Split the fasta file into individual files
awk '/^>/ {if(s) {print s > f; close(f)} f=sprintf("%s/%s.fasta", "'${split_fasta_dir}'", substr($0,2)); s=""} {s=s sprintf("%s\n",$0)} END {if(s) print s > f; close(f)}' "${original_fasta}"

# Function to run AlphaFold on a single split fasta file
run_alphafold() {
    local fasta_file=$1
    local gpu_id=$2 # GPU ID to use for this job

    cmd_args="--fasta_paths=${fasta_file}
    --output_dir=${project_root}outputs/$(basename ${fasta_file} .fasta)
    --db_preset=full_dbs
    --data_dir=${data_root}
    --max_template_date=2023-03-20
    --gpu_devices=${gpu_id}"

    CUDA_VISIBLE_DEVICES=${gpu_id} alphafold ${cmd_args}
}

export -f run_alphafold # Export the function for parallel or xargs

# Round-robin GPU assignment starts here
gpu_ids=(0 1 2 3) # Array of GPU IDs
num_gpus=${#gpu_ids[@]} # Number of GPUs

i=0
for fasta_file in ${split_fasta_dir}/*.fasta; do
    gpu_id=${gpu_ids[i % num_gpus]}
    run_alphafold "$fasta_file" "$gpu_id" &
    let i++
done
wait # Wait for all background jobs to finish

# module load 2022
# module load AlphaFold/2.3.1-foss-2022a-CUDA-11.7.0

# # export $PATH

# HHBLITS_PATH="/sw/arch/RHEL8/EB_production/2022/software/HH-suite/3.3.0-gompi-2022a/bin/hhblits"
# JACKHMMER_PATH="/sw/arch/RHEL8/EB_production/2022/software/HMMER/3.3.2-gompi-2022a/bin/jackhmmer"

# export PATH=$HHBLITS_DIR:$JACKHMMER_DIR:$PATH

# # Use the executables with their full paths

# $JACKHMMER_PATH -A ./project/alphafold/outputs/msa/uniref90_hits.sto --noali --F1 0.0005 --F2 5e-05 --F3 5e-07 --incE 0.0001 -E 0.0001 --cpu 8 -N 1 ./project/alphafold/inputs/6M0J.fasta /projects/2/managed_datasets/AlphaFold/uniref90/uniref90.fasta
# $JACKHMMER_PATH ./project/alphafold/outputs/msa/mgnify_output.sto --noali --F1 0.0005 --F2 5e-05 --F3 5e-07 --incE 0.0001 -E 0.0001 --cpu 8 -N 1 ./project/alphafold/inputs/6MOJ.fasta /projects/2/managed_datasets/AlphaFold/mgnify/mgy_clusters_2022_05.fa
# $HHBLITS_PATH -i ./project/alphafold/outputs/msa/6M0J.fasta -cpu 4 -oa3m ./project/alphafold/outputs/msa/bfd_uniref_hits.a3m -n 3 -e 0.001 -maxseq 1000000 -realign_max 100000 -maxfilt 100000 -min_prefilter_hits 1000 -d /projects/2/managed_datasets/AlphaFold/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt -d /projects/2/managed_datasets/AlphaFold/uniref30/UniRef30_2021_03

# source /etc/profile.d/modules.sh