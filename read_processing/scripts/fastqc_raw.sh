#!/bin/zsh

mkdir -p output/check_fastq/
fastqc input/*.fq.gz -o output/check_fastq/ -t 10 --nogroup

