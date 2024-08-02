#!/bin/bash

# Directory containing your PDB and PQR files
PDB_DIR="/home/max/stayahead/out_tmp/esmfold/controls"

# Directory to store output
OUTPUT_DIR="/home/max/stayahead/snellius2/outputs/apbs/test/esm_controls"

# Loop through all PDB files in the directory
for pdb_file in $PDB_DIR/*.pdb; do
    # Generate the PQR file from the PDB file using pdb2pqr
    pqr_file="${pdb_file%.pdb}.pqr"
    pdb2pqr --ff=AMBER --nodebump --noopt --apbs-input="${pdb_file%.pdb}.in" $pdb_file $pqr_file

    # Construct the input filename for APBS
    input_file="${pdb_file%.pdb}.in"

    # Run APBS and redirect output to a text file
    output_file="$OUTPUT_DIR/$(basename ${pdb_file%.pdb}).txt"
    apbs $input_file > $output_file
done

echo "APBS processing complete for all files."
