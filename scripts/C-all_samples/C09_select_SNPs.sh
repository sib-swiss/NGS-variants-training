#!/usr/bin/env bash

cd ~/workdir

gatk SelectVariants \
--variant results/variants/trio.vcf \
--select-type-to-include SNP \
--output results/variants/trio.SNP.vcf
