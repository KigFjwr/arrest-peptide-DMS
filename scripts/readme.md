# NGS Read Processing Scripts for DMS

This repository provides scripts for processing Illumina NGS reads used in Deep Mutational Scanning (DMS) experiments.
The workflow includes MD5 check, demultiplexing, read quality filtering/merging, pattern generation, fitness calculation, and plotting.

⸻

## Usage Overview

Step 0. Preparation  
- Place paired-end read files (read_1.fq.gz, read_2.fq.gz) and the MD5 checksum file in the input folder. 
- Add information to 'sample_sheet.csv'.   

⸻

## Step 1. MD5 Check

Script: check_md5sum.R  
Description: Verify integrity of downloaded read files.   
Dependencies: R (tidyverse)  

```bash
Rscript scripts/check_md5sum.R input/read_1.fq.gz input/read_2.fq.gz input/md5.txt
```

⸻

## Step 2. Demultiplexing

Script: Demultiplexing.sh  
Description: Split reads by sample using fastp   
Dependencies: fastp  
Note: The folder fastq_separated will be generated and can be deleted after use.   

```bash
zsh scripts/Demultiplexing.sh input/read_1.fq.gz input/read_2.fq.gz input/sample_sheet.csv
```

⸻

## Step 3. Quality Check & Read Merging

Script: MergeReads.sh    
Description: Assess read quality and merge paired-end reads.     
Dependencies: fastp    

```bash
zsh scripts/MergeReads.sh -s input/sample_sheet.csv -o <output_dir>
# with custom fastp options:
zsh scripts/MergeReads.sh -s input/sample_sheet.csv -o <output_dir> -q 20 -u 30 -e 25 -n 5
```

⸻

## Step 4. Generate DMS Pattern List

Script: PrepPattern_NNN.R or PrepPattern_NNK.R    
Description: Create a list of all possible DMS mutant sequences.     
Dependencies: R (tidyverse, Biostrings, gtools)    

```bash
Rscript scripts/PrepPattern_NNK.R <target_name> <target_sequence> <output_name>
```

⸻

## Step 5. Calculate Growth Rate

Script: calc_change.R    
Description: Compute growth rate change from merged read counts.     
Dependencies: R (tidyverse, ShortRead)    

```bash
Rscript scripts/calc_change.R <pattern_list> <5_seq> <3_seq> <merged_reads_dir> <DMS_name> <cutoff> <sample_sheet>
```

⸻

## Step 6. Plot Results

Script: summarize.R    
Description: Complute relative fitness and generate summary plots of DMS results.     
Dependencies: R (tidyverse, Biostrings, ggthemes, patchwork)    

```bash
Rscript scripts/summarize.R <DMS_name> <region_start> <region_end> <growth_rate_csv> <pattern_list> <output_suffix>
```

⸻
