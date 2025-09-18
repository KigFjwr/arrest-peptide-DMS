#!/usr/bin/zsh
set -e
set -u
set -o pipefail

# arguments
# argument $1: read 1 of paired-end reads
read_1=$1
# argument $2: read 2 of paired-end reads
read_2=$2

mkdir -p output/fastp/

fastp -i "$read_1" -I "$read_2" -h output/fastp/report_rawreads.html -j output/fastp/report_rawreads.json -w 8