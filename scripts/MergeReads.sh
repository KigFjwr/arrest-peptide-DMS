#!/usr/bin/zsh

set -e
set -u
set -o pipefail

#-------------------------------------------

# Arguments
# デフォルト値
sample_sheet="input/sample_sheet.csv"
output_dir="qc1"
qualified_quality_phred=20
unqualified_percent_limit=30
average_qual=25
n_base_limit=5


# getopts でオプションを処理
while getopts "s:o:q:u:e:n:" OPT; do
  case "$OPT" in
    s) sample_sheet="$OPTARG" ;;
    o) output_dir="$OPTARG" ;;
    q) qualified_quality_phred="$OPTARG" ;;
    u) unqualified_percent_limit="$OPTARG" ;;
    e) average_qual="$OPTARG" ;;
    n) n_base_limit="$OPTARG" ;;
  esac
done


shift $((OPTIND - 1))  # 残りの位置引数を処理
# 確認用出力
echo "sample_sheet: $sample_sheet"
echo "output_dir: $output_dir"
echo "qualified_quality_phred: $qualified_quality_phred"
echo "unqualified_percent_limit: $unqualified_percent_limit"
echo "average_qual: $average_qual"
echo "n_base_limit: $n_base_limit"

#-------------------------------------------

# 0. Prepare directories
mkdir -p output/fastq_merged/$output_dir


# 1. read processing
while IFS=, read sample_number bc5 bc7 sample_name read_exp lib_len target|| [ -n "${read_exp}" ]; do 
      # If the first character is "#" then jump to the next loop, considering that line as a comment
     if [ ${sample_number:0:1} = "#" ]; then 
        continue
     fi
  
  echo output/fastq_demultiplexed/$sample_name/*1.fq.gz | read input_read_1
  echo output/fastq_demultiplexed/$sample_name/*2.fq.gz | read input_read_2
  echo output/fastq_merged/$output_dir/${sample_name}_read_1.fq.gz | read output_read_1
  echo output/fastq_merged/$output_dir/${sample_name}_read_2.fq.gz | read output_read_2
  echo output/fastq_merged/$output_dir/${sample_name}_merged.fq.gz | read output_read_m
  echo output/fastq_merged/$output_dir/${sample_name}_unpair_1.fq.gz | read output_unpair_1
  echo output/fastq_merged/$output_dir/${sample_name}_unpair_2.fq.gz | read output_unpair_2

  fastp \
  -m \
  -i $input_read_1 \
  -I $input_read_2 \
  -o $output_read_1 \
  -O $output_read_2 \
  --merged_out $output_read_m \
  --unpaired1 $output_unpair_1 \
  --unpaired2 $output_unpair_2 \
  -h output/fastq_merged/$output_dir/report_${sample_name}.html \
  -j output/fastq_merged/$output_dir/report_${sample_name}.json \
  -q $qualified_quality_phred \
  -u $unqualified_percent_limit \
  -e $average_qual \
  -n $n_base_limit \
  -5 \
  -3 \
  -W 4 \
  -M 20 \
  -l $lib_len \
  -w 8

done < $sample_sheet

# last update
# 2025/10/08
