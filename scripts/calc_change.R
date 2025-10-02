# last update: 2025/10/01

# input -------------------------------------------------------------------

args <- commandArgs(trailingOnly = T)

default_args <- c("", "", "", "", "20", '8', 'input/sample_sheet.csv')   # default setting
default_flg <- is.na(args[1:7])
args[default_flg] <- default_args[default_flg]  


# argument 1: input file name of DMS-pattern list
# example: input_pattern <- 'output/ptn_ApdA_S109toP123_NNK.csv'
input_pattern <- args[1]

# argument 2: input sequence (5'->3' direction) which locate upstream of DMS-target region
# example: input_region_5f <- 'CCCGCGCTCCTGGGTGTGCACGGTCACGGTGATCTGCCTCTCTTCGGCACCGTCCCGCACGGACCGGCGTCCACGACGCTGCTCCGCCTCAGCGAGCTCCACGACGAAGCAGCTCCCGCC' # 5'->3'
input_region_5f <- args[2]

# argument 3: input sequence (5'->3' direction) which locate downstream of DMS-target region
# example: input_region_3f <- 'TTCCCGGTAGCCGTCGACTATAAAGACGACGACGACAAA' # 5'->3'
input_region_3f <- args[3]

# argument 4: suffix of the output file name
# example: output_suffix <- 'ApdA_S109-P123_PY79_qc20'
output_suffix <- args[4]


# optional

# argument 5: phred quality score used
# default: qc <- '20'
qc <- args[5]

# argument 6: read cutoff score
# default: read_cutoff = 10
read_cutoff <- as.numeric(args[6])

# argument 7: sample_sheet file
# default: sample_sheet <- 'input/sample_sheet.csv'
sample_sheet <- args[7]


# load package ------------------------------------------------------------

library(tidyverse)
library(ShortRead)

# source ------------------------------------------------------------------

# load input --------------------------------------------------------------

# load sample sheet
df_sample_sheet <- read_csv(sample_sheet) %>% dplyr::rename(sample_number = `#sample_number`) %>%
  bind_cols(
    bc7_RevCom = as.character(
      Biostrings::DNAStringSet(.$bc7) %>% 
        Biostrings::reverseComplement()
    )
  )

# load parttern sheet
df_pattern <- read_csv(input_pattern) %>% 
  bind_cols(
    seq_r = as.character(
      Biostrings::DNAStringSet(.$seq) %>% 
        Biostrings::reverseComplement()
    )
  )


# input read files list
read_files <- fs::dir_ls(path = 'output/fastp_qc/', recurse = T) %>% 
  str_subset(qc) %>% # select files that contain reads filterd by phred score indicated
  str_subset('merged') %>% # select pair-end-merged file
  print()


# load fastq file
fq_read <- readFastq(read_files)
df_read <- tibble(
  seq_raw = ShortRead::sread(fq_read) %>% as.character(),
  name = ShortRead::id(fq_read) %>% as.character() %>% str_extract('[[:graph:]]+')
)


# prepare output file names
output_rpm <- str_c('output/calc/CalcRPM_', output_suffix, '.csv')
output_fc <- str_c('output/calc/CalcGR_', output_suffix, '.csv')
output_filter_out <- str_c('output/calc/CalcGR_', output_suffix, '_removed.csv')
output_missing <- str_c('output/calc/CalcGR_', output_suffix, '_missing.csv')


# prepare output directory
if(!dir.exists('output/calc/')){
  dir.create('output/calc/')
}



# calculate RPM -----------------------------------------------------------
# RPM: Reads per million 


df_read_processed <- df_read %>% 
  
  # filtering

  ## Select reads which has perfect region_5'
  dplyr::filter(str_detect(seq_raw, input_region_5f)) %>%
  
  ## Select reads which has perfect region_3'
  dplyr::filter(str_detect(seq_raw, input_region_3f)) %>%

  
  # Sequence processing
  
  dplyr::mutate(
    
    # extract BC7 sequences
    bc7_RevCom = str_sub( 
      seq_raw,
      str_locate(seq_raw, input_region_3f)[,2] +1,
      str_locate(seq_raw, input_region_3f)[,2] +3
    ),
    
    # Extract DMS target region
    seq_target = str_sub( 
      seq_raw,
      str_locate(seq_raw, input_region_5f)[,2] +1,
      str_locate(seq_raw, input_region_3f)[1] -1
    )
  ) %>% 
  
  # combine sample information
  dplyr::left_join(
    df_sample_sheet %>% 
      dplyr::select(sample_number, bc7, bc7_RevCom, sample_name), by = 'bc7_RevCom'
  ) %>% 
  dplyr::select(sample_number, sample_name, bc7, name, everything()) %>%
  dplyr::filter(!is.na(bc7))


# mapping to df_pattern
df_read_mapped <- df_pattern %>%
  left_join(df_read_processed, by = c('seq' = 'seq_target')) %>% 
  dplyr::filter(!is.na(sample_number)) %>% 
  print()


# calculate size factor for RPM
df_SizeFactor <- df_read_mapped %>% 
  dplyr::count(sample_number) %>%
  dplyr::transmute(
    sample_number,
    # calculate size factors for RPM
    sf=1000000/n
  ) %>% 
  print()



# calculate RPM
df_tmp_RPM <- df_read_mapped %>% 
  # count reads per pattern
  dplyr::group_by(sample_number, pattern_number) %>% 
  dplyr::reframe(reads = n()) %>% 
  
  # calculate RPM
  dplyr::left_join(df_SizeFactor, by = join_by(sample_number)) %>% 
  dplyr::mutate(RPM = reads * sf) %>% 
  
  # Replace RPM with ‘NA’ if the read counts at t1 or t2 are below the cutoff
  dplyr::mutate(
    RPM = if_else(
      reads >= read_cutoff, 
      RPM,
      NA
    )
  )


# formatting
df_RPM <- df_pattern %>%
  left_join(
    df_sample_sheet,
    by = 'target',
    relationship = "many-to-many"
  ) %>% 
  dplyr::left_join(df_tmp_RPM, by = join_by(sample_number, pattern_number)) %>% 
  dplyr::select(sample_number, sample_name, organism, bc7, bc5, everything()) %>%
  dplyr::arrange(sample_number, desc(RPM)) %>% 
  ungroup() %>% 
  print()

df_RPM %>% 
  write_csv(output_rpm)



# calc change -------------------------------------------------------------




df_t1 <- df_RPM %>% 
  dplyr::filter(time == 't1') %>% 
  dplyr::select(
    pattern_number,
    rep_selection,
    t1_reads = reads,
    t1_cfu = cfu, 
    t1_culture_min = culture_min, 
    t1_RPM = RPM,
  )


df_change <- df_RPM %>% 
  dplyr::filter(time == 't2') %>% 
  dplyr::rename(
    t2_reads = reads,
    t2_cfu = cfu, 
    t2_culture_min = culture_min, 
    t2_RPM = RPM
  ) %>% 
  dplyr::left_join(df_t1, by = join_by(pattern_number, rep_selection)) %>% 
  group_by(sample_number) %>% 
  dplyr::mutate(
    relative_Fitness = log2(t2_RPM /t1_RPM ),
    rFitness_per_hour = log2(t2_RPM /t1_RPM ) / (t2_culture_min/60),
    Growth_Rate = log2(t2_cfu/t1_cfu * t2_RPM /t1_RPM ) / (t2_culture_min/60),
  ) %>% 
  ungroup() %>% 
  print()


df_change %>% 
  dplyr::select(target, organism, pattern_number, wtmt, mt_res, codon, seq, seq_aa, HasStopCodon, antibiotic, conc_antibiotic, rep_cell_lib, rep_selection, rep_DNA_lib,  t1_reads, t2_reads, t1_cfu, t2_cfu, t2_culture_min, t1_RPM, t2_RPM, relative_Fitness, rFitness_per_hour, Growth_Rate) %>%
  dplyr::arrange(pattern_number, rep_selection, conc_antibiotic) %>%
  write_csv(output_fc)



# optional information ----------------------------------------------------

# output removed raws
df_change %>% 
  filter(is.na(Growth_Rate)) %>% 
  dplyr::select(target, organism, pattern_number, wtmt, mt_res, codon, seq, seq_aa, HasStopCodon, antibiotic, conc_antibiotic, rep_cell_lib, rep_selection, rep_DNA_lib, t1_cfu, t2_cfu, t2_culture_min, t1_reads, t2_reads, t1_RPM, t2_RPM, relative_Fitness, rFitness_per_hour, Growth_Rate) %>% 
  write_csv(output_filter_out)


# output missing codon/aa information
df_codon_pattern <- df_pattern %>% 
  group_by(seq_aa) %>% 
  dplyr::count(name = 'n_codon_pattern')

df_change %>% 
  dplyr::filter(is.na(Growth_Rate)) %>% 
  dplyr::select(sample_number, sample_name, seq_aa, codon) %>% 
  dplyr::distinct() %>% 
  add_count(sample_number, sample_name, seq_aa, name = 'n_missing_codon') %>% 
  left_join(df_codon_pattern) %>% 
  select(-codon) %>% 
  dplyr::distinct() %>% 
  mutate(
    n_in_library = n_codon_pattern - n_missing_codon,
    InLibrary = case_when(
      n_in_library > 0 ~ 'missing_codon',
      n_in_library == 0 ~ 'missing_aa',
    )
  ) %>% 
  tidyr::complete(sample_name, seq_aa) %>% 
  select(sample_number, everything()) %>% 
  arrange(sample_number, seq_aa) %>% 
  write_csv(output_missing)





# last update: 2025/10/01
