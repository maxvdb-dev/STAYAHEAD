#!/bin/bash

# Navigate to the directory containing the zip files
cd test

# Find the zip file with the lowest number (assumes filenames are like ACE2_{number}.zip)
lowest_zip=$(ls ACE2_*.zip | sort -V | head -n 1)

# Extract the number from the filename
batch_number=$(echo $lowest_zip | sed -E 's/ACE2_([0-9]+)\.zip/\1/')

mkdir tmp$batch_number

cp $lowest_zip tmp$batch_number

rm $lowest_zip

# Create a new directory for this batch
mkdir batch$batch_number

# Unzip the file into the newly created directory
unzip tmp$batch_number/$lowest_zip -d batch$batch_number

# Delete the zipfile
rm -r tmp$batch_number

echo "Processed batch $batch_number."

