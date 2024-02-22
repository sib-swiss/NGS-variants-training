#!/usr/bin/env bash

cd ~/project

gatk VariantsToTable \
--variant results/variants/mother.HC.vcf \
--fields CHROM -F POS -F TYPE -GF GT \
--output results/variants/mother.HC.table
