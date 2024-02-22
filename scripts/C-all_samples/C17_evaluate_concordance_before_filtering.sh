#!/usr/bin/env bash

cd ~/project

gatk Concordance \
--evaluation results/variants/mother.trio.vcf \
--truth data/variants/NA12878.vcf.gz \
--intervals chr20:10018000-10220000 \
--summary results/variants/concordance.mother.trio
