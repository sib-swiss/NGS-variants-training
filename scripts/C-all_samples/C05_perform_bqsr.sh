#!/usr/bin/env bash

cd ~/workdir

for SAMPLE in mother father son
do
  gatk BaseRecalibrator \
  --reference data/reference/Homo_sapiens.GRCh38.dna.chromosome.20.fa \
  --input results/alignments/"$SAMPLE".rg.md.bam \
  --known-sites data/variants/GCF.38.filtered.renamed.vcf \
  --known-sites data/variants/1000g_gold_standard.indels.filtered.vcf \
  --output bqsr/"$SAMPLE".recal.table

  gatk ApplyBQSR \
  --input results/alignments/"$SAMPLE".rg.md.bam \
  --bqsr-recal-file bqsr/"$SAMPLE".recal.table \
  --output bqsr/"$SAMPLE".recal.bam
done
