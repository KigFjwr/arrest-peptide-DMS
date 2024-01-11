#!/bin/zsh

# input
score=20


cd /Users/kfsci1/Dropbox/BioInfo/proj_2023/p100_DMS_packaging/dms_analysis
mkdir -p output/fastp_qc/qc$score

## Loading the sample sheet
input_sample_sheet=input/sample_sheet.csv

while IFS=, read sample_number org bc5 bc7 sample_name read_exp || [ -n "${read_exp}" ]; do 
      # If the first character is "#" then jump to the next loop, considering that line as a comment
     if [ ${sample_number:0:1} = "#" ]; then 
        continue
     fi
  
  echo input/$sample_name/*1.fq.gz | read input_read_1
  echo input/$sample_name/*2.fq.gz | read input_read_2
  echo output/fastp_qc/qc$score/${sample_name}_read_1_qc${score}.fq.gz | read output_read_1
  echo output/fastp_qc/qc$score/${sample_name}_read_2_qc${score}.fq.gz | read output_read_2

  fastp -i $input_read_1 -I $input_read_2 \
  -o $output_read_1 -O $output_read_2 \
  -h output/fastp_qc/qc$score/report_${sample_name}.html \
  -j output/fastp_qc/qc$score/report_${sample_name}.json \
  -q $score \
  -e 35 \
  -l 20 \
  -w 8

done < $input_sample_sheet



