#!/usr/bin/env bash

cd ~/project/results

for SAMPLE in mother father son
do
    gatk MarkDuplicates \
    --INPUT alignments/"$SAMPLE".rg.bam \
    --OUTPUT alignments/"$SAMPLE".rg.md.bam \
    --METRICS_FILE alignments/marked_dup_metrics_"$SAMPLE".txt 
done
