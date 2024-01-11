#!/bin/zsh

# 0. Preparing the directory
mkdir -p output/fastq_sepalated/
mkdir -p output/fastq_sepalated/tmp/
mkdir -p output/fastq_demultiplexed/


# 1. Using "seqkit split" to split read_2 for each region with a barcode in the sequence
# argument $1: input file, read_2
seqkit split $1 -r 7:9 -O output/fastq_sepalated 



# 2. From the files split by barcode, rename only the necessary ones to the sample names listed in the pre-prepared sample sheet

## Loading the sample sheet
input_sample_sheet=input/sample_sheet_read_proc.csv

while IFS=, read sample_number bc5 bc7 sample_name read_exp || [ -n "${read_exp}" ]; do # Ensure the ability to read until the last line even if there is no newline code at the end of the final line

    # If the first character is "#" then jump to the next loop, considering that line as a comment
     if [ ${sample_number:0:1} = "#" ]; then 
        continue
     fi
     
     # Setting up the source file for copying
     echo output/fastq_sepalated/*$bc7.fq.gz | read copy_from
     # Setting up the destination file name for copying
     echo output/fastq_sepalated/tmp/${sample_name}.fq.gz | read copy_to
     # print
     echo $copy_from $copy_to
     # copy
     cp $copy_from $copy_to

done < $input_sample_sheet



# 3. Using "seqkit pair", prepare read_1 with its paired read_2
while IFS=, read sample_number bc5 bc7 sample_name read_exp || [ -n "${read_exp}" ]; do

    # If the first character is "#" then jump to the next loop, considering that line as a comment
     if [ ${sample_number:0:1} = "#" ]; then 
        continue
     fi

    mkdir -p output/fastq_demultiplexed/${sample_name}

    # Configuring the read_2 file
    echo output/fastq_sepalated/tmp/${sample_name}.fq.gz | read file_read_2
    # Configuring the read_1 file
    echo input/*_1.fq.gz | read file_read_1
    echo $file_read_1 $file_read_2
    
    # Matching the pairs of paired-end reads with "seqkit pair".
    seqkit pair -1 $file_read_1 -2 $file_read_2 -O output/fastq_demultiplexed/${sample_name}

    # Rename the files to distinguish between read_1 and read_2.
    mv output/fastq_demultiplexed/${sample_name}/${file_read_1##*/} output/fastq_demultiplexed/${sample_name}/${sample_name}_read_1.fq.gz
    mv output/fastq_demultiplexed/${sample_name}/${file_read_2##*/} output/fastq_demultiplexed/${sample_name}/${sample_name}_read_2.fq.gz


done < $input_sample_sheet
