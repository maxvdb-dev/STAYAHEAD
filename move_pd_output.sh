#!/bin/bash

# Define the base directory and the new output directory
NEW_OUTPUT_DIR="pd_output_ba1_100"

# Create the new output directory
# mkdir -p "$NEW_OUTPUT_DIR"

# List of output directories to process
OUTPUT_DIRS=("output_ba1_100")

# Loop through each output directory
for OUTPUT_DIR in "${OUTPUT_DIRS[@]}"; do
    FULL_OUTPUT_PATH="$OUTPUT_DIR"

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
                cp "$SUBDIR"/*ranked_0* "$NEW_SUBDIR"/
                
                # Check if ranked_debug.json exists and process it
                if [ -f "$SUBDIR/ranking_debug.json" ]; then
                    # Extract the best model name using jq
                    BEST_MODEL=$(jq -r '.order[0]' "$SUBDIR/ranking_debug.json")

                    # # Construct the corresponding pkl file name
                    PKL_FILE="$SUBDIR/result_${BEST_MODEL}.pkl"

                    # # Check if the pkl file exists and copy it
                    if [ -f "$PKL_FILE" ]; then
                        cp "$PKL_FILE" "$NEW_SUBDIR"/
                    else
                        echo "File $PKL_FILE does not exist."
                    fi

                    # Optionally copy the ranking_debug.json file if needed
                    cp "$SUBDIR/ranking_debug.json" "$NEW_SUBDIR"/
                else
                    echo "File $SUBDIR/ranking_debug.json does not exist."
                fi
            fi
        done
    else
        echo "Directory $FULL_OUTPUT_PATH does not exist."
    fi
done

echo "File extraction and copying complete."