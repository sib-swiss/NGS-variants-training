#!/usr/bin/env bash

cd ~/workdir/results/alignments
samtools flagstat mother.sam > mother.sam.flagstat
