#!/bin/bash
#SBATCH --job-name=MMSEQS2_job
#SBATCH --nodes=1
#SBATCH --time=06:00:00
#SBATCH --error=mmseqs2/logs/run_mmseqs_jobs_%A_%a_err.txt
#SBATCH --output=mmseqs2/logs/run_mmseqs_jobs_%A_%a_out.txt
#SBATCH --partition=genoa
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem 256G
#SBATCH --array=8-15%2

echo "Local scratch space path: $TMPDIR"

NODE_NAME=$(hostname)
echo "Node name: $NODE_NAME"

BATCH_NUM=$(printf "%d" $SLURM_ARRAY_TASK_ID)

TMP_OUTPUT="$TMPDIR/output${BATCH_NUM}"

mkdir -p $TMP_OUTPUT

# Ensure directories are created successfully
if [ ! -d "$TMP_OUTPUT" ]; then
    echo "Error creating TMP_OUTPUT directory: $TMP_OUTPUT" >&2
    exit 1
fi

ZIP_PATH="./mmseqs2/fastas/ACE2_100000/ACE2_${BATCH_NUM}.zip"
BATCH_DIR="$TMPDIR/batch${BATCH_NUM}"

mkdir -p $BATCH_DIR

# Ensure directories are created successfully
if [ ! -d "$BATCH_DIR" ]; then
    echo "Error creating BATCH_DIR directory: $BATCH_DIR" >&2
    exit 1
fi

unzip $ZIP_PATH -d $BATCH_DIR

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

ZIP_OUTPUT="$TMP_OUTPUT.zip"
HOME_ZIP="$HOME/stayahead/mmseqs2/outputs/2step/output${BATCH_NUM}.zip"

zip -r $ZIP_OUTPUT $TMP_OUTPUT

mv $ZIP_OUTPUT $HOME_ZIP


# echo "Zipping output directory..."
# zip -r $TMPDIR/output_${BATCH_NUM}.zip -C $TMP_OUTPUT .

# echo "Copying zipped output to final destination..."
# cp $TMPDIR/output_${BATCH_NUM}.zip $HOME/stayahead/mmseqs2/outputs/2step/

# rm -r $TMP_OUTPUT

# /scratch-shared/tmp.W3tCMjupp5