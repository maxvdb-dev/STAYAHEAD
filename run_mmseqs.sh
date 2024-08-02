#!/bin/bash
#SBATCH --job-name=MMSEQS2_job
#SBATCH --nodes=1
#SBATCH --time=05:00:00
#SBATCH --error=mmseqs2/logs/run_mmseqs_jobs_%A_%a_err.txt
#SBATCH --output=mmseqs2/logs/run_mmseqs_jobs_%A_%a_out.txt
#SBATCH --partition=genoa
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem 256G
#SBATCH --array=5-6%2

BATCH_NUM=$(printf "%d" $SLURM_ARRAY_TASK_ID)

TMP_DIR="/home/mvandenboom/scratch-shared/tmp.W3tCMjupp5"

TMP_OUTPUT="$TMP_DIR/output${BATCH_NUM}"

mkdir -p $TMP_OUTPUT

# Ensure directories are created successfully
if [ ! -d "$TMP_OUTPUT" ]; then
    echo "Error creating TMP_OUTPUT directory: $TMP_OUTPUT" >&2
    exit 1
fi

# ZIP_PATH="./mmseqs2/fastas/test/ACE2_${BATCH_NUM}.zip"
BATCH_DIR="$TMP_DIR/batch${BATCH_NUM}"

mkdir -p $BATCH_DIR

# Ensure directories are created successfully
if [ ! -d "$BATCH_DIR" ]; then
    echo "Error creating BATCH_DIR directory: $BATCH_DIR" >&2
    exit 1
fi

# SCRIPT_PATH="./mmseqs2/fastas/unzip_batch.sh"

# unzip $ZIP_PATH -d $BATCH_DIR

# cp -r ./mmseqs2/fastas/ACE2_100000/batch4 $TMPDIR/batch4
# echo "Copied batch to $TMPDIR/batch4"
# rm -r ./mmseqs2/fastas/ACE2_100000/batch4
# source colabfold env
source ./project/colabfold/venv/bin/activate

# load uniref indices into mem
vmtouch -t -l -d -w /projects/2/managed_datasets/AlphaFold_mmseqs2/uniref30_2302_db.idx

# Define command arguments
cmd_args="--use-env 0 
--use-templates 0 
--db-load-mode 2 
--mmseqs ./mmseqs2/mmseqs/bin/mmseqs 
--threads 64
--db1 uniref30_2302_db $BATCH_DIR
/projects/2/managed_datasets/AlphaFold_mmseqs2/
$TMP_OUTPUT"

# search, .a3m file will be stored in ./output_test
colabfold_search ${cmd_args}

echo "Zipping output directory..."
zip -r $TMP_DIR/output_${BATCH_NUM}.zip -C $TMP_OUTPUT .

echo "Copying zipped output to final destination..."
cp $TMP_DIR/output_${BATCH_NUM}.zip $HOME/stayahead/mmseqs2/outputs/2step/

rm -r $TMP_OUTPUT

# /scratch-shared/tmp.W3tCMjupp5