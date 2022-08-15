#!/usr/bin/env bash

cd ~/workdir/data

gatk IndexFeatureFile --input variants/1000g_gold_standard.indels.filtered.vcf
gatk IndexFeatureFile --input variants/GCF.38.filtered.renamed.vcf
