# NGS Read Processing and Analysis Scripts for Deep Mutational Scanning (DMS)

This repository contains scripts used for NGS data processing and analysis
in a deep mutational scanning (DMS) study.

The analysis starts from demultiplexed FASTQ files deposited in the
DNA Data Bank of Japan (DRA).
Scripts and parameters are provided to enable reproducible re-analysis
from the public data.

---

## Data source and starting point

Sequencing reads deposited in DRA are already demultiplexed and correspond
to individual experimental conditions.
These demultiplexed FASTQ files were used as the starting point of all
analyses described in this repository.

---

## Directory structure (recommended)
```
.
├── input/
│   ├── *_read_1.fq.gz
│   ├── *_read_2.fq.gz
│   ├── plot_order.csv
│   └── sample_sheet.csv
├── scripts/
│   ├── MergeReads.sh
│   ├── PrepPattern_NNK.R
│   ├── PrepPattern_NNN.R
│   ├── calc_change.R
│   └── summarize.R
└── output/
    ├── fastq_demultiplexed/
    ├── ptn/
    ├── calc/
    └── figures/
```

---

## 1. Preparation

- Download demultiplexed FASTQ files from DRA and place them in `output/fastq_demultiplexed`
- Prepare `sample_sheet.csv` and place it in `input/`

Note:
The last column of `sample_sheet.csv` must be named `read_exp`.

---

## 2. Read quality filtering and paired-end merging

Script:
- `MergeReads.sh`

Description:
Reads are quality-filtered and merged using fastp.
Quality thresholds are explicitly defined to ensure reproducibility.

Dependency:
- fastp

Example:
```
cd /path/to/project
zsh scripts/MergeReads.sh \
  -s input/sample_sheet.csv \
  -o qc1 \
  -q 20 \
  -u 30 \
  -e 25 \
  -n 5
```

⸻

3. Preparation of DMS mutant pattern list

Scripts:
	•	PrepPattern_NNK.R
	•	PrepPattern_NNN.R

Description:
Generate a list of expected mutant sequences for DMS libraries
based on a wild-type nucleotide sequence.

Note:
Depending on the library design, either NNK or NNN is used.
These scripts are alternatives and should not be used sequentially
for the same dataset.

Dependencies:
	•	R
	•	tidyverse
	•	Biostrings
	•	gtools

Example:
```
cd /path/to/project
Rscript scripts/PrepPattern_NNK.R \
  TARGET_NAME \
  WILD_TYPE_SEQUENCE \
  output/ptn/ptn_TARGET_NNK.csv
```

⸻

4. Fitness calculation

Script:
	•	calc_change.R

Description:
Calculates fitness changes for each mutant based on read counts
after quality filtering and merging.

Example:
```
cd /path/to/project
Rscript scripts/calc_change.R \
  output/ptn/ptn_TARGET_NNK.csv \
  UPSTREAM_SEQUENCE \
  DOWNSTREAM_SEQUENCE \
  qc1 \
  TARGET_NAME \
  8 \
  input/sample_sheet.csv
```

⸻

5. Visualization and summary

Script:
	•	summarize.R

Description:
Summarizes DMS results and generates plots.

Example:
```
cd /path/to/project
Rscript scripts/summarize.R \
  TARGET_NAME \
  START_POSITION \
  END_POSITION \
  output/calc/qc1/GR_TARGET_NAME.csv \
  output/ptn/ptn_TARGET_NNK.csv \
  qc1
```

⸻

Software requirements (tested versions)
	•	fastp v0.24.0
	•	SeqKit v2.10.0
	•	R v4.5.0
	•	R packages:
	•	tidyverse v2.0.0
	•	Biostrings v2.76.0
	•	ShortRead v1.66.0
	•	gtools v3.9.5
	•	ggthemes v5.1.0
	•	patchwork v1.3.2

⸻

Versioning

The version of the code used in the published manuscript is tagged on GitHub.
Later commits may include refactoring or additional features not used
in the original analysis.

---
