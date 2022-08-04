#!/usr/bin/env bash

cd ~/workdir/results

for SAMPLE in mother father son
do
    samtools index alignments/"$SAMPLE".rg.md.bam
done  
