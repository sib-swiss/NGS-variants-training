#!/usr/bin/env bash

cd ~/workdir/results

gatk MergeVcfs \
--INPUT variants/trio.SNP.filtered.vcf \
--INPUT variants/trio.INDEL.filtered.vcf \
--OUTPUT variants/trio.filtered.vcf
