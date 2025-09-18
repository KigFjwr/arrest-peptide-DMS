# NGS read processing scripts for DMS

## 0. preparation
- put read files (read_1 and read_2) into the input folder
- put the md5 file supplied from NGS service to the input folder
- prepare sample_sheet_read_proc.csv and put it into input folder  
note: last column (named read_exp) is required in sample_sheet_read_proc.csv  


## 1. check MD5
script: check_md5sum.R  
description: a code to check downloaded files  
required: input/[read_1].fq.gz input/[read_2].fq.gz, input/md5.txt  

```usage example
$ Rscript scripts/check_md5sum.R [input/read1.fq.gz] [input/read2.fq.gz] [input/md5.txt]
```

## 2. check read quality
script: fastqc_raw.sh  
description: a code to check read quality  

```
$ cd /path/to/read_processing
$ zsh scripts/fastqc_raw.sh
```

## 3. demultiplexing
script: demultiplexing_seqkit.sh  
description: a code for demultiplexing reads  
required: input/sample_sheet_read_proc.csv  
note: fastq_sepalated folder will be created, but there is no need to keep it for a long time. Please delete it as needed.  

```
$ zsh scripts/demultiplexing_seqkit.sh [input/read_1.fq.gz] [input/read_2.fq.gz]
```


## 4. Prep patten list
script: PrepPattern_NNN.R or PrepPattern_NNK.R
description: a code for  
required: DMS target sequence  

```
Rscript PrepPattern_NNN.R <seq_name> <seq> <output>

```

---
## required (checekd version)
- fastp v0.23.4
- SeqKit v2.5.1
- R version v4.3.0
- R packages
  - tidyverse, v2.0.0
  - Biostrings, v2.72.1
  - ShorRead, v1.62.0
