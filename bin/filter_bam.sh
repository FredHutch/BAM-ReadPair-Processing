#!/bin/bash

set -e

IN_BAM=$1
OUT_BAM=$2
THREADS=$3

samtools view -h -b -q 30 -F 0x4 -f 0x3 -@ $THREADS $IN_BAM > ${OUT_BAM}
samtools index -@ $THREADS ${OUT_BAM}
