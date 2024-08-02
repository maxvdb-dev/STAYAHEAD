#!/bin/bash

# Navigate to the batch32 directory
cd ba1_100

# Loop through each subdirectory
for dir in */
do
    # Check if the msas/msas directory exists
    if [ -d "$dir$dir/msas" ]; then
        # Move all files from the msas/msas directory two levels up
        mv "$dir$dir/msas/"* "$dir$dir/"
    fi
done

# Return to the original directory
cd ..