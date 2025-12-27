#!/usr/bin/zsh
set -e
set -u
set -o pipefail
setopt extended_glob null_glob

#-------------------------------------------

# Arguments

# argument $1: input file, read_1
file_read_1=$1
# argument $2: input file, read_2
file_read_2=$2
# argument $3: input file, sample sheet
sample_sheet=${3:-input/sample_sheet.csv}


#-------------------------------------------


# 0. Prepare directories
mkdir -p output/fastq_separated/
mkdir -p output/fastq_separated/tmp/
mkdir -p output/fastq_demultiplexed/


# 1. Using "seqkit split" to split read_2 for each region with a barcode in the sequence
seqkit split $file_read_2 -r 7:9 -O output/fastq_separated 



# 2. From the files split by barcode, rename only the necessary ones to the sample names listed in the pre-prepared sample sheet

while IFS=, read sample_number bc5 bc7 sample_name read_exp || [ -n "${read_exp}" ]; do # Ensure the ability to read until the last line even if there is no newline code at the end of the final line

    # If the first character is "#" then jump to the next loop, considering that line as a comment
     if [ ${sample_number:0:1} = "#" ]; then 
        continue
     fi
     
    # find source file: accept .fq.gz or .fastq.gz
    matches=(output/fastq_separated/*${bc7}.(fq|fastq).gz(N))

    if (( ${#matches} == 0 )); then
      echo "ERROR: No file matched bc7=${bc7} in output/fastq_separated/" >&2
      exit 1
    elif (( ${#matches} > 1 )); then
      echo "ERROR: Multiple files matched bc7=${bc7}:" >&2
      printf '%s\n' "${matches[@]}" >&2
      exit 1
    fi

    copy_from="${matches[1]}"
    copy_to="output/fastq_separated/tmp/${sample_name}.fq.gz"

     # print
     echo $copy_from $copy_to
     # copy
     cp $copy_from $copy_to

done < $sample_sheet



# 3. Using "seqkit pair", prepare read_1 with its paired read_2
while IFS=, read sample_number bc5 bc7 sample_name read_exp || [ -n "${read_exp}" ]; do

    # If the first character is "#" then jump to the next loop, considering that line as a comment
     if [ ${sample_number:0:1} = "#" ]; then 
        continue
     fi

    mkdir -p output/fastq_demultiplexed/${sample_name}

    # Configuring the read_2 file
    echo output/fastq_separated/tmp/${sample_name}.fq.gz | read file_read_2
    echo $file_read_1 $file_read_2
    
    # Matching the pairs of paired-end reads with "seqkit pair".
    seqkit pair -1 $file_read_1 -2 $file_read_2 -O output/fastq_demultiplexed/${sample_name}

    # Rename the files to distinguish between read_1 and read_2.
    mv output/fastq_demultiplexed/${sample_name}/${file_read_1##*/} output/fastq_demultiplexed/${sample_name}/${sample_name}_read_1.fq.gz
    mv output/fastq_demultiplexed/${sample_name}/${file_read_2##*/} output/fastq_demultiplexed/${sample_name}/${sample_name}_read_2.fq.gz


done < $sample_sheet

## version 2512227