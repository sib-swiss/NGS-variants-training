#!/usr/bin/env bash

cd ~/project

for SAMPLE in mother father son
do
    bwa mem data/reference/Homo_sapiens.GRCh38.dna.chromosome.20.fa \
    data/fastq/"$SAMPLE"_R1.fastq.gz \
    data/fastq/"$SAMPLE"_R2.fastq.gz \
    | samtools sort \
    | samtools view -bh > results/alignments/"$SAMPLE".bam
done
