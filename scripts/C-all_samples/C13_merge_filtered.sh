#!/usr/bin/env bash

cd ~/workdir/results/variants

gatk MergeVcfs \
--INPUT trio.SNP.filtered.vcf \
--INPUT trio.INDEL.filtered.vcf \
--OUTPUT trio.filtered.vcf