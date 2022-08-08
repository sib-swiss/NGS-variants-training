#!/usr/bin/env bash

cd ~/workdir/results

gatk SelectVariants \
--variant variants/trio.filtered.vcf \
--sample-name mother \
--exclude-non-variants \
--remove-unused-alternates \
--output variants/mother.trio.filtered.vcf
