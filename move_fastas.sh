#!/bin/bash

# Set the source and destination root directories
src_root="1step100"
dest_root="1step100_fastas"

mkdir -p "$dest_root"

# Loop through each batch directory in the source root directory
for batch_dir in "$src_root"/batch*/; do
    # Ensure the directory exists and is not empty
    if [ -d "$batch_dir" ]; then
        echo "Processing batch directory: $batch_dir"
        # Copy all .fasta files from the project directory to the destination directory
        cp "$batch_dir"/*.fasta "$dest_root"/
    fi
done

echo "All .fasta files have been copied to $dest_root"
