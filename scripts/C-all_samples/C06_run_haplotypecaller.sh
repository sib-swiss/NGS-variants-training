#!/usr/bin/env bash

cd ~/project

for SAMPLE in mother father son
do
  gatk HaplotypeCaller \
  --reference data/reference/Homo_sapiens.GRCh38.dna.chromosome.20.fa \
  --input results/bqsr/"$SAMPLE".recal.bam \
  --output results/variants/"$SAMPLE".HC.g.vcf \
  --bam-output results/variants/"$SAMPLE".phased.bam \
  --intervals chr20:10018000-10220000 \
  --emit-ref-confidence GVCF
done
