#!/bin/bash

# Starting directory
root_dir="."

# Find all the 'msas' directories under each 'batch' directory
find "${root_dir}" -type d -name "msas" | while read -r msas_dir; do
    # Get the full path of the parent directory of 'msas'
    parent_dir=$(dirname "${msas_dir}")
    
    # Extract the name of the parent directory
    parent_name=$(basename "${parent_dir}")
    
    # Construct the new directory path, which is a subdirectory of the parent with the same name
    new_dir="${parent_dir}/${parent_name}"
    
    # Create the new directory
    mkdir -p "${new_dir}"
    
    # Move the 'msas' directory into the newly created directory
    mv "${msas_dir}" "${new_dir}/"
    
    echo "Moved ${msas_dir} to ${new_dir}/msas"
done
