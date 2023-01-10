#!/usr/bin/env bash

cd ~/workdir/results/variants

gatk VariantFiltration \
--variant trio.SNP.vcf \
--filter-expression "QD < 2.0"              --filter-name "QD2" \
--filter-expression "QUAL < 30.0"           --filter-name "QUAL30" \
--filter-expression "SOR > 3.0"             --filter-name "SOR3" \
--filter-expression "FS > 60.0"             --filter-name "FS60" \
--filter-expression "MQ < 40.0"             --filter-name "MQ40" \
--filter-expression "MQRankSum < -12.5"     --filter-name "MQRankSum-12.5" \
--filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
--output trio.SNP.filtered.vcf
