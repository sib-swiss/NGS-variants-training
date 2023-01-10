#!/usr/bin/env bash

cd ~/workdir/results/variants

gatk SelectVariants \
--variant trio.vcf \
--sample-name mother \
--exclude-non-variants \
--remove-unused-alternates \
--select-type-to-include INDEL \
--select-type-to-include SNP \
--output mother.trio.vcf
