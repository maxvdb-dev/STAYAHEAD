#!/bin/bash
#SBATCH --job-name=compute_metrics
#SBATCH --nodes=1
#SBATCH --time=01:00:00 
#SBATCH --partition=genoa
#SBATCH --tasks-per-node=1 
#SBATCH --cpus-per-task=18 

# module load Python/3.9.5-GCCcore-10.3.0
module load 2022
module load Python/3.10.4-GCCcore-11.3.0
module load Biopython/1.79-foss-2022a

python3.10 -c "from Bio import PDB; print('Biopython is correctly installed!')"

root_dir="./project/alphafold/outputs/msa/1step/batch11_15"
ref_file="./project/alphafold/inputs/6m0j.pdb"
out_dir="./metrics/alphafold/1step50/batch11_15/"
# root_dir=/project/alphafold/outputs/ds4
# ref_file=/project/alphafold/inputs/6m0j.pdb

# Loop through each subdirectory in the specified directory
for subdir in "$root_dir"/*; do
    # Construct the path to the ranked_0.pdb file within the current subdirectory
    pdb_file="$subdir/$(basename "$subdir")/ranked_0.pdb"
    
    # Check if the pdb file exists before calling the Python script
    if [ -f "$pdb_file" ]; then
        # Call your Python script here, passing the PDB file path as an argument
        # Replace 'your_script.py' with the actual name of your Python script
        python3.10 comp_analysis.py -i "$pdb_file" -o "$out_dir" -r "$ref_file"
    else
        echo "File not found: $pdb_file"
    fi
done