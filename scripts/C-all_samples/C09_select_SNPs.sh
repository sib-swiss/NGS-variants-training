#!/usr/bin/env bash

cd ~/project

gatk SelectVariants \
--variant results/variants/trio.vcf \
--select-type-to-include SNP \
--output results/variants/trio.SNP.vcf
