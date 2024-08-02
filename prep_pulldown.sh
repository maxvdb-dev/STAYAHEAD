#!/bin/bash

# Set the source and destination root directories
src_root="zips"
dest_root="pd100_missing"
a_dir="ACE2-ectodomain"

candidates_file="${dest_root}/candidates_100_missing.txt"

# Create the destination root directory if it doesn't exist
mkdir -p "$dest_root"

# Initialize or clear the candidates.txt file
> "$candidates_file"

# Loop through each batch directory
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

                # Define the path to the msas directory
                msas_dir="${subdir}msas"

                # Check if the msas directory exists
                if [ -d "$msas_dir" ]; then
                    # Extract the project directory name
                    project_name=$(basename "$project_dir")
                    # Define the destination directory path
                    dest_dir="${dest_root}/${project_name}"
                    dest_dir="${dest_root}"
                    features_file="${subdir}features.pkl"
                    # Create the destination directory
                    mkdir -p "$dest_dir"

                    # Copy the msas directory to the destination directory and rename it
                    cp -r "$msas_dir" "$dest_dir/${project_name}"
                    cp -r "$a_dir" "$dest_dir"
                    cp "$features_file" "$dest_dir/${project_name}.pkl"
                    echo "$project_name" >> "$candidates_file"
                fi
            done
        done
    fi
done

echo "All msas directories have been processed."
