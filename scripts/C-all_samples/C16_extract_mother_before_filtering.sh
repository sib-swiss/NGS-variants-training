#!/usr/bin/env bash

cd ~/workdir/results

gatk SelectVariants \
--variant variants/trio.vcf \
--sample-name mother \
--exclude-non-variants \
--remove-unused-alternates \
--select-type-to-include INDEL \
--select-type-to-include SNP \
--output variants/mother.trio.vcf