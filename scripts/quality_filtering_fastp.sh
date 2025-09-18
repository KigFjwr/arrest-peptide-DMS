#!/bin/zsh

# input
score=20


mkdir -p output/fastp_qc/qc$score

## Load the sample sheet
# argument $1: sample_sheet_library.csv
input_sample_sheet=$1

while IFS=, read sample_number org object bc5 bc7 sample_name read_exp || [ -n "${read_exp}" ]; do 
      # If the first character is "#" then jump to the next loop, considering that line as a comment
     if [ ${sample_number:0:1} = "#" ]; then 
        continue
     fi
  
  echo output/fastq_demultiplexed/$sample_name/*1.fq.gz | read input_read_1
  echo output/fastq_demultiplexed/$sample_name/*2.fq.gz | read input_read_2
  echo output/fastp_qc/qc$score/${sample_name}_read_1_qc${score}.fq.gz | read output_read_1
  echo output/fastp_qc/qc$score/${sample_name}_read_2_qc${score}.fq.gz | read output_read_2
  echo output/fastp_qc/qc$score/${sample_name}_merged_qc${score}.fq.gz | read output_read_m
  echo output/fastp_qc/qc$score/${sample_name}_unpair_1_qc${score}.fq.gz | read output_unpair_1
  echo output/fastp_qc/qc$score/${sample_name}_unpair_2_qc${score}.fq.gz | read output_unpair_2

  fastp \
  -q $score \
  -i $input_read_1 \
  -I $input_read_2 \
  -o $output_read_1 \
  -O $output_read_2 \
  --merged_out $output_read_m \
  --unpaired1 $output_unpair_1 \
  --unpaired2 $output_unpair_2 \
  -h output/fastp_qc/qc$score/report_${sample_name}.html \
  -j output/fastp_qc/qc$score/report_${sample_name}.json \
  -c \
  -m \
  -e 35 \
  -l 20 \
  -w 8

done < $input_sample_sheet

