#!/usr/bin/env nextflow

// Enable DSL 2 syntax
nextflow.enable.dsl = 2

process filter_bam {
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
        tuple path(bam), path(bai)
    
    output:
        path "output/*.bam", emit: bam
        path "output/*_readpair_counts.csv", emit: csv

    """#!/bin/bash
set -e

# Create the output directory
mkdir -p output

# Run the bam_count.py script
bam_count.py "${bam}" "output/${bam}"
    """
}

process merge_csv {
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
        path "output/"
    
    output:
        path "output/merged_reference_count.csv"

    """#!/bin/bash
set -e

# Run the merge_csv.py script
merge_csv.py output
    """
}

workflow {

    if ( params.help ){
        log.info"""
        BAM-ReadPair-Processing

        Required Arguments:
            --indir         Folder containing BAM files to process
            --outdir        Folder where output BAM files will be written

        """
        exit 0
    }

    // Make a channel with all of the BAM files
    // and then run the filter_bam process on them
    Channel
        .fromPath("${params.indir}/*.bam")
        .map { it -> [it, file("${it}.bai", checkIfExists: true)]}
        .ifEmpty { error "No files found matching the pattern ${params.indir}/*.bam" }
        | filter_bam

    // Combine all of the CSVs
    merge_csv(
        filter_bam.out.csv.toSortedList()
    )
}
