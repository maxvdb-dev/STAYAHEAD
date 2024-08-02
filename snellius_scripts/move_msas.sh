#!/bin/bash

# Starting directory
root_dir="ba1_100"

# Find directories containing the MSA files
find "${root_dir}" -type d | while read -r dir_path; do
    # Count the number of files in the directory
    num_files=$(find "${dir_path}" -maxdepth 1 -type f | wc -l)
    
    # Proceed if the directory contains files
    if [ "${num_files}" -gt 0 ]; then
        # Path for the new 'msas' directory
        msas_dir="${dir_path}/msas"
        
        # Create the 'msas' directory if it doesn't exist
        mkdir -p "${msas_dir}"
        
        # Move all files in the directory to 'msas'
        # Using find to select only files avoids moving the newly created 'msas' directory
        find "${dir_path}" -maxdepth 1 -type f -exec mv {} "${msas_dir}/" \;
        
        echo "Moved ${num_files} files to ${msas_dir}"
    fi
done
