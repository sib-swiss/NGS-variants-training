#!/usr/bin/env bash

cd ~/workdir/results

gatk Concordance \
--evaluation variants/mother.trio.filtered.vcf \
--truth data/variants/NA12878.vcf.gz \
--intervals chr20:10018000-10220000 \
--summary variants/concordance.mother.trio.filtered
