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


def skip_alignment(alignment, max_distance=1000) -> bool:
    """Return True if the alignment should be skipped"""

    if not alignment.is_proper_pair:
        return True
    elif alignment.is_read2:
        return True
    elif alignment.next_reference_id != alignment.reference_id:
        return True
    elif calc_distance(alignment) > max_distance:
        return True
    else:
        return False


def calc_distance(alignment) -> int:
    """Calculate the distance between reads in a pair."""

    return abs(alignment.next_reference_start - alignment.reference_start)


def process_and_write_bam(input_bam_file, output_bam_file):

    # input bam file
    # (note: using the 'with' syntax ensures proper file handle closure)
    with pysam.AlignmentFile(input_bam_file, "rb", threads=30) as seqbam:

        # use defaultdict to init readpair_dict
        readpair_dict = defaultdict(list)

        # skip if not qualified
        for alignment in seqbam:

            # Check if the alignment should be skipped
            if skip_alignment(alignment):
                continue

            # The 'alignment' should always be read1
            assert alignment.is_read1

            # mate will be used to store r2 alignment
            mate = seqbam.mate(alignment)
            assert mate.is_read2

            readpair_dict[alignment.query_name].append({
                "AS": int(alignment.get_tag('AS')) + int(mate.get_tag('AS')),
                "reference_name": alignment.reference_name,
                "alignment": alignment,
                "mate": mate
            })

        # Copy the header to use in the output
        bam_header = seqbam.header

    # output bam file
    with pysam.AlignmentFile(output_bam_file, 'wb', header=bam_header) as outfile:

        # init a dict to count the READ PAIRS per reference
        readpair_count = defaultdict(int)

        for data in readpair_dict.values():
            # Get the highest alignment score
            max_score_data = max(data, key=lambda x: x['AS'])

            # count the read pair if there is only one with the highest AS
            if len(data) == 1:
                readpair_count[max_score_data['reference_name']] += 1

                # Write out the alignments to the output BAM file
                outfile.write(max_score_data["alignment"])
                outfile.write(max_score_data["mate"])

    return readpair_count


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_bam_file', type=str, help='input bam file')
    parser.add_argument('output_bam_file', type=str, help='output bam file')

    args = parser.parse_args()

    sample_name = os.path.basename(args.output_bam_file)[:-len('.bam')]
    output_path = os.path.dirname(args.output_bam_file)

    readpair_counts = process_and_write_bam(
        args.input_bam_file,
        args.output_bam_file
    )

    readpair_counts_df = pd.DataFrame(readpair_counts.items(), columns=["reference_name", sample_name])
    # calculate the total number of mapped pairs based on other reference count
    total_mapped_pairs = readpair_counts_df[sample_name].sum()

    # create a header df with 'Total_mapped_pairs'
    total_mapped_df = pd.DataFrame({"reference_name": ["Total_mapped_pairs"], 
                                    sample_name: [total_mapped_pairs]})

    # then put the header and reference counts together
    readpair_counts_df = pd.concat([total_mapped_df, readpair_counts_df], ignore_index=True)

    readpair_counts_df.to_csv(
        output_path + '/' + sample_name + "_readpair_counts.csv", index=False
    )


if __name__ == "__main__":
    main()
