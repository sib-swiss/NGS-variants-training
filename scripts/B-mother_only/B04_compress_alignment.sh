#!/usr/bin/env bash

cd ~/workdir/results

samtools view -bh alignments/mother.sorted.sam > alignments/mother.bam
