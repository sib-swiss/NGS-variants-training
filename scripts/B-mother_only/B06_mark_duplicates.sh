#!/usr/bin/env bash

cd ~/project/results

gatk MarkDuplicates \
--INPUT alignments/mother.rg.bam \
--OUTPUT alignments/mother.rg.md.bam \
--METRICS_FILE alignments/marked_dup_metrics_mother.txt 
