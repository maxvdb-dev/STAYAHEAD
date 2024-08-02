#!/bin/bash

# Starting directory
root_dir="ba1_100"

# Find all 'msas' directories
find "${root_dir}" -type d -name "msas" | while read -r msas_dir; do
    # Parent directory of the 'msas' directory
    parent_dir=$(dirname "${msas_dir}")
    
    # Count the number of files in the 'msas' directory
    num_files=$(find "${msas_dir}" -maxdepth 1 -type f | wc -l)
    
    # Proceed if the 'msas' directory contains files
    if [ "${num_files}" -gt 0 ]; then
        # Move all files from 'msas' directory to the parent directory
        find "${msas_dir}" -maxdepth 1 -type f -exec mv {} "${parent_dir}/" \;
        
        echo "Moved ${num_files} files from ${msas_dir} to ${parent_dir}"
    fi
    
    # Remove the empty 'msas' directory
    rmdir "${msas_dir}"
    echo "Removed directory ${msas_dir}"
done
