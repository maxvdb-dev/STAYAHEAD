#!/bin/bash

# Set the source and destination root directories
src_root="all_batches_100"
dest_root="1step_all_100"

mkdir -p "$dest_root"

for batch_dir in "$src_root"/batch*/; do
    # Ensure the directory exists and is not empty
    if [ -d "$batch_dir" ]; then
        echo "Processing batch directory: $batch_dir"
        # Loop through each project directory in the batch directory
        for project_dir in "$batch_dir"/*/; do
            # Loop through the duplicated subdirectory inside each project directory
            for subdir in "$project_dir"/*/; do
                # Print the subdirectory being processed
                echo "Currently processing subdir: $subdir"

                project_name=$(basename "$project_dir")

                cp "${subdir}ranked_0.pdb" "$dest_root/${project_name}.pdb"
            done
        done
    fi
done