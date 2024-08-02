#!/bin/bash

# Define the base directory and the new output directory
BASE_DIR="/scratch-shared/tmp.KJVsRs6W2S"
NEW_OUTPUT_DIR="$BASE_DIR/af_output_controls"

# Create the new output directory
mkdir -p "$NEW_OUTPUT_DIR"

# List of output directories to process
OUTPUT_DIRS=("af_controls")

# Loop through each output directory
for OUTPUT_DIR in "${OUTPUT_DIRS[@]}"; do
    FULL_OUTPUT_PATH="$BASE_DIR/$OUTPUT_DIR"

    # Check if the directory exists
    if [ -d "$FULL_OUTPUT_PATH" ]; then
        # Loop through each subdirectory in the output directory
        for SUBDIR in "$FULL_OUTPUT_PATH"/*; do
            if [ -d "$SUBDIR" ]; then
                # Extract the directory name
                DIR_NAME=$(basename "$SUBDIR")

                # Create the corresponding directory in the new output directory
                NEW_SUBDIR="$NEW_OUTPUT_DIR/$DIR_NAME"
                mkdir -p "$NEW_SUBDIR"

                # Copy the required files to the new directory
                cp "$SUBDIR"/"$DIR_NAME"/*ranked_0* "$NEW_SUBDIR"/
            fi
        done
    else
        echo "Directory $FULL_OUTPUT_PATH does not exist."
    fi
done

echo "File extraction and copying complete."