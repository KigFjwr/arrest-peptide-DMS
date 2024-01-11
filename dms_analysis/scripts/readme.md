# DMS analysis

## 0. preparation
- demultiplexed read files , in input folder
- prepare sample_sheet.csv
-   
note:   


## 1. Prepare pattern list
script: PrepPattern_NNN.R etc  
description: a code to prepare a list of all mutation pattern  
required:   
note: ランダム（N）などの入れ方によって、使用するスクリプトを選択する。
- NNNの場合はPrepPattern_NNN.R
- NNXの場合は
```usage example
$ cd /path/to/dms_analysis
$ Rscript scripts/PrepPattern_NNN.R [peptide_name] [seq] [output file name]
```

## 2. Read filteration
script:   
required: fastp  
description: a code to check read quality  
note: fastp is run by x64 architecture. 
```
$ cd /path/to/dms_analysis
$ zsh scripts/quality_filtering_fastp.sh
```

```
# fastp command
$ fastp -i $input_read_1 -I $input_read_2 \
  -o $output_read_1 -O $output_read_2 \
  -h output/fastp_qc/qc$score/report_${sample_name}.html \
  -j output/fastp_qc/qc$score/report_${sample_name}.json \
  -q $score \
  -e 35 \
  -l 20 \
  -w 8
```

### fastp options;  
input/output
- -i, --in1: input, read_1  
- -I, --in2: input, read_2
- -o, --out1: output, read_1
- -O, --out2: output, read_2

Quality
- -q, --qualified_quality_phred: the quality value that a base is qualified. Default 15 means phred quality >=Q15 is qualified. (int [=15])
- -e, --average_qual: if one read's average quality score <avg_qual, then this read/pair is discarded. Default 0 means no requirement (int [=0])

Trimming
- -n, --n_base_limit: if one read's number of N base is >n_base_limit, then this read/pair is discarded. Default is 5 (int [=5])
- -l, --length_required: reads shorter than length_required will be discarded, default is 15. (int [=15])

Others
- -w: worker thread number, default is 3 (int [=3])

## 3. Calculate Fitness
script: calc_fitness.R  
description:  
required:  
note:   

```
$ cd /path/to/dms_analysis
$ Rscript scripts/calc_change.R [DMS-pattern list] [seq_5'] [seq_3'] [target] [q] [co] [sample_sheet]
```

arguments
- [DMS-pattern list]
- [seq_5'] = sequence (5'->3' direction) which locate upstream of DMS-target region (e.g. ATTGCGCCAGAGATGGCTCCTTCCCTGCCGGTGGCGAAAACCAGAATTGCGCGTCTCCCATCCTGT)
- [seq_3'] = sequence (5'->3' direction) which locate downstream of DMS-target region (e.g. GCGGCGGGAGCCTTCCTTGACTATAAAGACGACGACGACAAA)
- [target] = target discreption, used as a suffix of the output file name (e.g. ApdP_Q126-P134_JM109)
- [q] = phred quality score used (default = 20)
- [co] = read cutoff score (default = 8)
- [sample_sheet] = sample_sheet.csv (default = input/sample_sheet.csv)



## 4. Plotting
script: plotting.R  
description:  
required:  
note:   
```
$ cd /path/to/dms_analysis
$ Rscript scripts/plotting.R [target] [Rstart] [Rend] [fitness list] [DMS-pattern list]
```

arguments
- [target] = target discreption, used as a suffix of the output file name (e.g. ApdP_Q126-P134_JM109)
- [Rstart] = start residue number of DMS target region (eg 126)
- [Rend] = last residue number of DMS target region (eg 134)
- [fitness list] = an (eg output/calc/CalcFC_[target].csv)
- [DMS-pattern list] = pattern list (eg ptn_[target]_NNN.csv)



---
## required (checekd version)
- R version v4.3.0
- R packages
  - tidyverse
  - Biostrings
  - gtools
- fastp v0.23.4
