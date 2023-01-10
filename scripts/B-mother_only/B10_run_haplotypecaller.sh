#!/usr/bin/env bash

cd ~/workdir
mkdir -p results/variants

gatk HaplotypeCaller \
--reference data/reference/Homo_sapiens.GRCh38.dna.chromosome.20.fa \
--input results/bqsr/mother.recal.bam \
--output results/variants/mother.HC.vcf \
--intervals chr20:10018000-10220000
