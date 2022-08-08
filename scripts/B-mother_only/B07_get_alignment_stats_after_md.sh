#!/usr/bin/env bash

cd ~/workdir/results/alignments

samtools flagstat mother.rg.md.bam > mother.rg.md.bam.flagstat
