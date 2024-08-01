#!/bin/bash

# Starting directory
root_dir="."

# Use find to locate all 'mgnify_output.sto' files under the specified directory structure
find "${root_dir}" -type f -name "mgnify_output.sto" | while read -r file_path; do
    # Directory of the current file
    dir=$(dirname "${file_path}")
    
    # New file name with 'mgnify_hits.sto'
    new_file_path="${dir}/mgnify_hits.sto"
    
    # Rename the file
    mv "${file_path}" "${new_file_path}"
    
    echo "Renamed ${file_path} to ${new_file_path}"
done
