#!/usr/bin/env bash

cd ~/project/results

gatk AddOrReplaceReadGroups \
--INPUT alignments/mother.bam \
--OUTPUT alignments/mother.rg.bam \
--RGLB lib1 \
--RGPU H0164.2.ALXX140820 \
--RGPL ILLUMINA \
--RGSM mother \
--RGID H0164.2
