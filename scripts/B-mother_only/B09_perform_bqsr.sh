#!/usr/bin/env bash

cd ~/workdir

mkdir -p results/bqsr

gatk BaseRecalibrator \
--reference data/reference/Homo_sapiens.GRCh38.dna.chromosome.20.fa \
--input results/alignments/mother.rg.md.bam \
--known-sites data/variants/GCF.38.filtered.renamed.vcf \
--known-sites data/variants/1000g_gold_standard.indels.filtered.vcf \
--output results/bqsr/mother.recal.table

gatk ApplyBQSR \
--input results/alignments/mother.rg.md.bam \
--bqsr-recal-file results/bqsr/mother.recal.table \
--output results/bqsr/mother.recal.bam
