#!/usr/bin/env bash

cd ~/workdir

gatk VariantFiltration \
--variant results/variants/trio.INDEL.vcf \
--filter-expression "QD < 2.0"                  --filter-name "QD2" \
--filter-expression "QUAL < 30.0"               --filter-name "QUAL30" \
--filter-expression "FS > 200.0"                --filter-name "FS200" \
--filter-expression "ReadPosRankSum < -20.0"    --filter-name "ReadPosRankSum-20" \
--output results/variants/trio.INDEL.filtered.vcf
