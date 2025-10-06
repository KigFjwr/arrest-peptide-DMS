# NGS read processing scripts for DMS

## 0. prep
- put read files (read_1 and read_2) into the input folder
- put the md5 file supplied from NGS service to the input folder
- prepare sample_sheet_read_proc.csv and put it into input folder  
note: last column (named read_exp) is required in sample_sheet_read_proc.csv  


## 1. check MD5
script: check_md5sum.R  
description: a code to check downloaded files  
required: input/[read_1].fq.gz, input/[read_2].fq.gz, input/md5.txt  
dependency: R, tidyverse  

```usage example
$ cd '/Users/kfsci4/Microb_ Physiol_Lab Dropbox/Fujiwara Keigo/BioInfo/proj_2025/p152_DMS/SecM_bs'
$ Rscript scripts/check_md5sum.R input/PY79_SecM3Y09B_1.fq.gz input/PY79_SecM3Y09B_2.fq.gz input/md5.txt
```


## 2. demultiplexing
script: demultiplexing_seqkit.sh.sh  
description: a code for demultiplexing reads  
required: input/sample_sheet_read_proc.csv  
note: fastq_sepalated folder will be created, but there is no need to keep it for a long time. Please delete it as needed.  

```
$ cd '/Users/kfsci4/Microb_ Physiol_Lab Dropbox/Fujiwara Keigo/BioInfo/proj_2025/p152_DMS/SecM_bs'
$ zsh scripts/Demultiplexing.sh input/PY79_SecM3Y09B_1.fq.gz input/PY79_SecM3Y09B_2.fq.gz input/sample_sheet.csv
```


## 3. check read quality and merge paired end reads
script: MergeReads.sh  
description: a code to check read quality  
required: input/[read_1].fq.gz, input/[read_2].fq.gz 
dependency: fastp  
```
$ cd '/Users/kfsci4/Microb_ Physiol_Lab Dropbox/Fujiwara Keigo/BioInfo/proj_2025/p152_DMS/SecM_bs'
$ zsh scripts/MergeReads.sh input/sample_sheet.csv 20
$ zsh scripts/MergeReads.sh input/sample_sheet.csv 30
```

## 4. prepare DMS pattten list
script: PrepPattern_NNN.R or PrepPattern_NNK.R  
description: codes to prepare DMS mutant sequence list  
required: DMS target sequence (wild-type)  
dependency: R, tidyverse, Biostrings, gtools

```
$ cd '/Users/kfsci4/Microb_ Physiol_Lab Dropbox/Fujiwara Keigo/BioInfo/proj_2025/p152_DMS/SecM_bs'  
# example, NNK: $ Rscript scripts/PrepPattern_NNK.R [target name] [target sequence] [output name]  
# example, NNN: $ Rscript scripts/PrepPattern_NNN.R [target name] [target sequence] [output name]  
$ Rscript scripts/PrepPattern_NNN.R SecM CCGTCTGAAAAGGGTTATCGCATTGATTATGCGCATTTTACCCCACAAGCAAAATTCAGCACGCCCGTCTGGATAAGCCAGGCGCAAGGCATCCGTGCTGGCCCT output/ptn_SecM_P132toP166_NNN.csv
```

## 5. calc fitness
script: calc_change.R  
description: codes to calculate fitness  
required:    
dependency: 
```
$ cd '/Users/kfsci4/Microb_ Physiol_Lab Dropbox/Fujiwara Keigo/BioInfo/proj_2025/p152_DMS/SecM_bs'  
$ Rscript scripts/calc_change.R output/ptn_SecM_P132toP166_NNN.csv CCGCAAACACTGCCCGTTGCTGAAGAATCTTTGCCTCTTCAGGCGCAACATCTTGCATTACTGGATACGCTCAGCGCGCTGCTGACCCAGGAAGGCACG CAACGCCTCACCGACTATAAAGACGACGACGACAAA SecM_P132toP166_PY79_qc20_co8 20 8 input/sample_sheet.csv
$ Rscript scripts/calc_change.R output/ptn_SecM_P132toP166_NNN.csv CCGCAAACACTGCCCGTTGCTGAAGAATCTTTGCCTCTTCAGGCGCAACATCTTGCATTACTGGATACGCTCAGCGCGCTGCTGACCCAGGAAGGCACG CAACGCCTCACCGACTATAAAGACGACGACGACAAA SecM_P132toP166_PY79_qc30_co8 30 8 input/sample_sheet.csv
$ Rscript scripts/calc_change.R output/ptn_SecM_P132toP166_NNN.csv CCGCAAACACTGCCCGTTGCTGAAGAATCTTTGCCTCTTCAGGCGCAACATCTTGCATTACTGGATACGCTCAGCGCGCTGCTGACCCAGGAAGGCACG CAACGCCTCACCGACTATAAAGACGACGACGACAAA SecM_P132toP166_PY79_qc20_co10 20 10 input/sample_sheet.csv
$ Rscript scripts/calc_change.R output/ptn_SecM_P132toP166_NNN.csv CCGCAAACACTGCCCGTTGCTGAAGAATCTTTGCCTCTTCAGGCGCAACATCTTGCATTACTGGATACGCTCAGCGCGCTGCTGACCCAGGAAGGCACG CAACGCCTCACCGACTATAAAGACGACGACGACAAA SecM_P132toP166_PY79_qc30_co10 30 10 input/sample_sheet.csv
```

## 6. plotting
script: summarize.R  
description: codes for plotting  
required:         
dependency:       
```
$ cd '/Users/kfsci4/Microb_ Physiol_Lab Dropbox/Fujiwara Keigo/BioInfo/proj_2025/p152_DMS/SecM_bs'  
$ Rscript scripts/summarize.R SecM_P132toP166_PY79 132 166 output/calc/CalcGR_SecM_P132toP166_PY79_qc20_co8.csv output/ptn_SecM_P132toP166_NNN.csv qc20_co8
$ Rscript scripts/summarize.R SecM_P132toP166_PY79 132 166 output/calc/CalcGR_SecM_P132toP166_PY79_qc30_co8.csv output/ptn_SecM_P132toP166_NNN.csv qc30_co8
$ Rscript scripts/summarize.R SecM_P132toP166_PY79 132 166 output/calc/CalcGR_SecM_P132toP166_PY79_qc20_co10.csv output/ptn_SecM_P132toP166_NNN.csv qc20_co10
$ Rscript scripts/summarize.R SecM_P132toP166_PY79 132 166 output/calc/CalcGR_SecM_P132toP166_PY79_qc30_co10.csv output/ptn_SecM_P132toP166_NNN.csv qc30_co10
```

---
## required (checekd version)
- R version v4.3.0
- R packages
  - tidyverse
  - Biostrings
  - gtools
- FastQC v0.12.1
- SeqKit v2.5.1
- 


