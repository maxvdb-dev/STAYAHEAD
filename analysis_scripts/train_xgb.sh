#!/bin/bash
#SBATCH --job-name=xgboost
#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --partition gpu
#SBATCH --gpus 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 18

#SBATCH -e logs/train_xgb_job_AF_bind_expr_err.txt
#SBATCH -o logs/train_xgb_job_AF_bind_expr_out.txt

## -e logs/train_xgb_job_ESM_bind_expr_deltas_err.txt
## -o logs/train_xgb_job_ESM_bind_expr_deltas_out.txt

module load 2023
module load Python/3.11.3-GCCcore-12.3.0
module load scikit-learn/1.3.1-gfbf-2023a
module load matplotlib/3.7.2-gfbf-2023a

python xgb_optuna.py -i ./inputs/AF/no_complex -o ./outputs/AF_train/1/bind_expr -p ./inputs/AF/no_complex/params.json -tt bind_expr -at alphafold
# python xgb_optuna.py -i ./inputs/AF/no_complex -o ./outputs/AF_train/1/bind_expr_deltas -p ./inputs/AF/no_complex/params.json -tt bind_expr_deltas -at alphafold
# python xgb_optuna.py -i ./inputs/ESM/no_complex -o ./outputs/train/ESM/bind_expr_deltas -p ./inputs/ESM/no_complex/params.json -tt bind_expr_deltas -at esmfold
# python xgb_optuna.py -i ./inputs/ESM/no_complex -o ./outputs/train/ESM/bind_expr -p ./inputs/ESM/no_complex/params.json -tt bind_expr -at esmfold