#!/usr/bin/env bash

cd ~/project/results

for SAMPLE in mother father son
do
    samtools index alignments/"$SAMPLE".rg.md.bam
done  
