#!/usr/bin/env python3
'''
This script will take one BWA-MEM bam file as the input, 
and output one bam file with only closely and properly aligned read pairs. 
The output also includes one read count csv and one reference count csv. 
Later, csv files from all input bam files will be integrated using Pandas.

USAGE:

python bam_count.py input_bam_file output_bam_file
'''

import os
import pysam
import pandas as pd
from collections import defaultdict
import argparse


def count_bam(input_bam_file, threads):
    # input bam file
    with pysam.AlignmentFile(
        input_bam_file,
        "rb",
        threads=int(threads)
    ) as seqbam:

        # init a dict to count the READ PAIRS per reference
        readpair_count = defaultdict(int)

        # No need to check for skip_alignment or max AS score due to proprocessing with Samtools
        for alignment in seqbam:
            # Read1 + Read2 = 1
            if alignment.is_read1:
                ref_name = alignment.reference_name + '_R1'
                readpair_count[ref_name] += 1
            elif alignment.is_read2:
                ref_name = alignment.reference_name + '_R2'
                readpair_count[ref_name] += 1

    return readpair_count


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_bam_file', type=str, help='input bam file')
    parser.add_argument('threads', type=str, help='input bam file')

    args = parser.parse_args()

    sample_name = os.path.basename(args.input_bam_file)[:-len('.bam')]

    readpair_counts = count_bam(args.input_bam_file, args.threads)

    readpair_counts_df = pd.DataFrame(
        readpair_counts.items(),
        columns=["reference_name", sample_name]
    )

    # calculate the total number of mapped pairs based on other reference count
    total_mapped_pairs = readpair_counts_df[sample_name].sum()

    # create a header df with 'Total_mapped_pairs'
    total_mapped_df = pd.DataFrame({"reference_name": ["Total_mapped_reads"], 
                                    sample_name: [total_mapped_pairs]})

    # then put the header and reference counts together
    readpair_counts_df = pd.concat(
        [total_mapped_df, readpair_counts_df],
        ignore_index=True
    )

    readpair_counts_df.to_csv(
        "output/" + sample_name + "_readpair_counts.csv",
        index=False
    )


if __name__ == "__main__":
    main()
