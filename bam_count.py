'''
This script will take one BWA-MEM bam file as the input, 
and output one bam file with only closely and properly aligned read pairs. 
The output also includes one read count csv and one reference count csv. 
Later, csv files from all input bam files will be integrated using Pandas.

USAGE:
python bam_count.py input_bam_file output_bam_file
'''

import pysam
import pandas as pd
from collections import defaultdict
import argparse

def process_and_write_bam(input_bam_file, output_bam_file):
    # input bam file
    seqbam = pysam.AlignmentFile(input_bam_file, "rb",threads=30)

    # use defaultdict to init readpair_dict
    readpair_dict = defaultdict(list)

    # skip if not qualified
    for alignment in seqbam:
        if not alignment.is_proper_pair or alignment.is_read2 or alignment.next_reference_id != alignment.reference_id or abs(alignment.next_reference_start - alignment.reference_start) > 1000:
            continue
        
        # mate will be used to store r2 alignment
        mate = seqbam.mate(alignment)
        readpair_dict[alignment.query_name].append({
            "AS": int(alignment.get_tag('AS')) + int(mate.get_tag('AS')),
            "reference_name": alignment.reference_name,
            "alignment": alignment,
            "mate": mate
        })

    # output bam file
    outfile = pysam.AlignmentFile(output_bam_file, 'wb', header=seqbam.header)

    # init a dict to count the READ PAIRS per reference
    readpair_count = defaultdict(int)

    #read1_count and read2_count will be also used to validate the result
    read1_count = 0
    read2_count = 0

    for readpair, data in readpair_dict.items():
        # sort the data and take the one with the highest AS
        max_score_data = sorted(data, key=lambda x: x['AS'], reverse=True)[0]

        # count the read pair if there is only one data with the highest AS
        if data.count(max_score_data) == 1:
            readpair_count[max_score_data['reference_name']] += 1
            outfile.write(max_score_data["alignment"])
            outfile.write(max_score_data["mate"])

            if max_score_data["alignment"].is_read1:
                read1_count += 1
            if max_score_data["mate"].is_read2:
                read2_count += 1

    seqbam.close()
    outfile.close()

    return readpair_count, read1_count, read2_count

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_bam_file', type=str, help='input bam file')
    parser.add_argument('output_bam_file', type=str, help='output bam file')

    sample_name = os.path.basename(args.output_bam_file).rstrip('.bam')
    output_path = os.path.dirname(args.output_bam_file)

    args = parser.parse_args()

    readpair_counts, read1_count, read2_count = process_and_write_bam(args.input_bam_file, args.output_bam_file)
    
    readpair_counts_df = pd.Series(readpair_counts)

    readpair_counts_df.to_csv(output_path + '/' + sample_name + "_readpair_counts.csv")

    read_counts_output_csv = output_path + '/' + sample_name + "_read_counts.csv"
    read_counts_output = open(read_counts_output_csv,'w')

    read_counts_output_csv.write('bam_file,read1_count,read2_count\n')
    read_counts_output_csv.write(os.path.basename(args.input_bam_file) + ',' + str(read1_count) + ',' + str(read2_count) + '\n')

if __name__ == "__main__":
    main()
