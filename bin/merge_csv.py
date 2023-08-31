#!/usr/bin/env python3
import pandas as pd
import os
import argparse


def merge_csv_files(directory_path):
    # csv_files is a list of count csv files for each sample
    csv_files = [
        f
        for f in os.listdir(directory_path)
        if f.endswith('readpair_counts.csv')
    ]

    # n will be used to count csv
    n = 0
    for csv_file in csv_files:
        # Read the current CSV into a DataFrame
        current_df = pd.read_csv(os.path.join(directory_path, csv_file))

        # init the merged_df
        if n == 0:
            merged_df = pd.DataFrame()
            merged_df = current_df
        else:
            # start from the second dataframe, merge
            merged_df = pd.merge(
                merged_df,
                current_df,
                on='reference_name',
                how='outer'
            )
        n += 1
    
    output_file = os.path.join(directory_path, 'merged_reference_count.csv')
    # Repalce NaN with zero
    merged_df = merged_df.fillna(0)
    
    # Reindex based on colnames
    # merged_df = merged_df.sort_index(axis=1)

    merged_df.to_csv(output_file, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Merge all read pair count csv file in the folder"
    )
    parser.add_argument(
        'directory_path',
        type=str,
        help='Path to the folder with read pair count csv files'
    )
    args = parser.parse_args()
    merge_csv_files(args.directory_path)

# usage: python merge_csv.py READ_PAIR_CSV_FOLDER
