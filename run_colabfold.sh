#!/bin/bash
#SBATCH --job-name=colabfold_job
#SBATCH --nodes=1
#SBATCH --time=05:00:00
#SBATCH --partition=gpu
#SBATCH --gpus=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

# Activate environment and load necessary modules
module load 2022
module load CUDA/11.7.0
module load SciPy-bundle/2022.05-foss-2022a
module load NCCL/2.12.12-GCCcore-11.3.0-CUDA-11.7.0
module load protobuf/3.19.4-GCCcore-11.3.0
module load protobuf-python/3.19.4-GCCcore-11.3.0
module load TensorFlow/2.11.0-foss-2022a-CUDA-11.7.0
module load Biopython/1.79-foss-2022a
module load HH-suite/3.3.0-gompi-2022a
module load HMMER/3.3.2-gompi-2022a
module load Kalign/3.3.5-GCCcore-11.3.0
module load jax/0.3.25-foss-2022a-CUDA-11.7.0
module load OpenMM/7.7.0-foss-2022a-CUDA-11.7.0 

source ./venv/bin/activate

# Define command arguments
cmd_args="./fastas/alt_variants_batch4
./output_alt_var/batch4/
--data /projects/2/managed_datasets/AlphaFold/ 
--num-models 3 
--num-recycle 1"

colabfold_batch ${cmd_args}