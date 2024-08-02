import argparse
import csv
import pandas as pd
import os
from pathlib import Path

# python /home/max/stayahead/snellius2/scripts/analysis/link_ab_escape.py -i 


def create_parser():
    parser = argparse.ArgumentParser(description="Compute structural metrics.")
    parser.add_argument(
        "-i", "--input_dir", help="Path to input directory", type=Path, required=True,
    )
    parser.add_argument(
        "-o", "--out_dir", help="Path to save directory for sequences", type=Path, required=True,
    )
    parser.add_argument(
        "-l", "--link_file", help="Path to data to link with input", type=Path, required=True,
    )
    return parser

def process_directory(input_dir, link_file, out_dir):
    link_df = pd.read_csv(link_file)
    for file in os.listdir(input_dir):
        if file.endswith(".csv"):
            df = pd.read_csv(input_dir / file)
            merged_df = pd.merge(link_df, df, left_on='sequence name', right_on='mut_position', how='inner')
            merged_df.to_csv(out_dir / file, index=False)

def main():
    parser = create_parser()
    args = parser.parse_args()
    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)
    process_directory(args.input_dir, args.out_dir, args.link_file)

if __name__ == "__main__":
    main()