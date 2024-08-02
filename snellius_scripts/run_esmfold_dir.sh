#!/bin/bash
#SBATCH --job-name=esmfold_job
#SBATCH --nodes=1
#SBATCH --time=20:00:00
#SBATCH --partition=gpu
#SBATCH --gpus=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18

# Activate environment and load necessary modules
source ./venv/bin/activate
module load 2022
module load cuDNN/8.6.0.163-CUDA-11.8.0
cd ./esmfold/

# Directory base paths
input_base="./inputs/1step"
output_base="./outputs/1step"

# Loop over directories in input_base
for input_dir in ${input_base}/*; do
    # Check if it is a directory
    if [ -d "${input_dir}" ]; then
        # Extract the batch name
        batch_name=$(basename ${input_dir})

        # Define the output directory for the current batch
        output_dir="${output_base}/${batch_name}/"

        # Create the output directory if it does not exist
        mkdir -p ${output_dir}

        # Define command arguments
        cmd_args="--fastas_folder ${input_dir} \
        --output_folder ${output_dir}"

        # Run ESMFold
        python esmfold.py ${cmd_args}
    fi
done
