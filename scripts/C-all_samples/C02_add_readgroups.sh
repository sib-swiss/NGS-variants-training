#!/usr/bin/env bash

cd ~/project/results

cat sample_rg_fields.txt | while read SAMPLE LB PU ID
do
    gatk AddOrReplaceReadGroups \
    --INPUT alignments/"$SAMPLE".bam \
    --OUTPUT alignments/"$SAMPLE".rg.bam \
    --RGLB "$LB" \
    --RGPU "$PU" \
    --RGPL ILLUMINA \
    --RGSM "$SAMPLE" \
    --RGID "$ID"
done 
