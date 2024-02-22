#!/usr/bin/env bash

cd ~/project/results

samtools view -bh alignments/mother.sorted.sam > alignments/mother.bam
