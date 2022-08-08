#!/usr/bin/env bash

cd ~/workdir

gatk GenomicsDBImport \
--variant results/variants/mother.HC.g.vcf \
--variant results/variants/father.HC.g.vcf \
--variant results/variants/son.HC.g.vcf \
--intervals chr20:10018000-10220000 \
--genomicsdb-workspace-path results/genomicsdb
