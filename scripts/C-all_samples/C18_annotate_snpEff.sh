#!/usr/bin/env bash

cd ~/workdir/results/variants

snpEff -Xmx4g \
-v \
GRCh38.99 \
trio.filtered.vcf \
> trio.filtered.snpeff.vcf
