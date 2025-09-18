
# input -------------------------------------------------------------------

args <- commandArgs(trailingOnly = T)

default_args <- c("", "", "", "", "20", '8', 'input/sample_sheet.csv')   # default setting
default_flg <- is.na(args[1:7])
args[default_flg] <- default_args[default_flg]  


# argument 1: input file name of DMS-pattern list
# example: input_pattern <- 'output/ptn_CdCliM_K38toK77_NNN.csv'
input_pattern <- args[1]

# argument 2: input sequence (5'->3' direction) which locate upstream of DMS-target region
# example: input_region_5f <- 'AGAGACCACATGGTCCTTCTTGAGTTTGTAACAGCTGCTGGGATTACACATGGCATGGATGAACTATACAAAAAAGACCTCTTAAATCATAAAATT' # ApdP(104-125) 5'->3'
input_region_5f <- args[2]

# argument 3: input sequence (5'->3' direction) which locate downstream of DMS-target region
# example: input_region_3f <- 'GACTATAAAGACGACGACGACAAA' # ApdP(135-140)-FLAG, 
input_region_3f <- args[3]

# argument 4: suffix of the output file name
# example: output_suffix <- 'CliM_K38-K77_PY79'
output_suffix <- args[4]


# optional

# argument 5: phred quality score used
# default: qc = '20'
qc <- args[5]

# argument 6: read cutoff score
# default: read_cutoff = 8
read_cutoff <- as.numeric(args[6])

# argument 7: sample_sheet file
# sample_sheet <- 'input/sample_sheet_library.csv'
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
  str_subset(qc) %>% # select files filterd by phred score indicated
  str_subset('merged') %>% # select pair-end-merged file
  print()


# load fastq file
fq_readR <- readFastq(read_files)
# only use read_2 data
df_readR <- tibble(
  seq_raw = ShortRead::sread(fq_readR) %>% as.character(),
  name = ShortRead::id(fq_readR) %>% as.character() %>% str_extract('[[:graph:]]+')
)


# prepare reverse complement sequence of 5'/3' region
input_region_5r = as.character(Biostrings::DNAStringSet(input_region_5f) %>% Biostrings::reverseComplement())
input_region_3r = as.character(Biostrings::DNAStringSet(input_region_3f) %>% Biostrings::reverseComplement())



# prepare output file names
output_rpm <- str_c('output/calc/CalcRPM_', output_suffix, '.csv')
output_fc <- str_c('output/calc/CalcFC_', output_suffix, '.csv')
output_filter_out <- str_c('output/calc/CalcFC_', output_suffix, '_cutoff.csv')
output_lacked <- str_c('output/calc/CalcFC_', output_suffix, '_lacked.csv')

# prepare output directory
if(!dir.exists('output/calc/')){
  dir.create('output/calc/')
}

# calculate RPM -----------------------------------------------------------
# RPM: Reads per million 


df_readR_processed <- df_readR %>% 
  
  # filtering
  # filter by length
  dplyr::filter(str_length(seq_raw) > 240) %>%

  ## Select reads which has perfect region_5'
  dplyr::filter(str_detect(seq_raw, input_region_5f)) %>%

  ## Select reads which has perfect region_3'
  dplyr::filter(str_detect(seq_raw, input_region_3f)) %>%

  # filtering
  ## filter by quality
  ## filter out reads in which 'N' found within (except for last residue)
  filter(
    !str_detect(
      seq_raw %>%
        str_sub(10, -2),
      'N'
    )
  ) %>%
  
  
  
  # Sequence processing
  
  dplyr::mutate(
    #bc7 = str_sub(seq_raw, 7, 9), 
    bc7_RevCom = str_sub( # extract BC7 sequences
      seq_raw,
      str_locate(seq_raw, input_region_3f)[,2] +1,
      str_locate(seq_raw, input_region_3f)[,2] +3
    ),
    #seq_gene = str_sub(seq_raw, 10, -1), # trim adapter region
    seq_gene = str_sub(# trim adapter region
      seq_raw, 
      str_locate(seq_raw, input_region_5f)[1],
      str_locate(seq_raw, input_region_3f)[,2]
    ), 
    seq_target = str_sub( # Extract scanning target region (ApdP(126-134))
      seq_raw,
      str_locate(seq_raw, input_region_5f)[,2] +1,
      str_locate(seq_raw, input_region_3f)[1] -1
    )
  ) %>% 
  
  # formatting
  dplyr::left_join(df_sample_sheet %>% dplyr::select(sample_number, bc7, bc7_RevCom, sample_name), by = 'bc7_RevCom') %>% 
  dplyr::select(sample_number, sample_name, bc7, name, everything()) %>%
  dplyr::filter(!is.na(bc7)) %>% 
  print()



# mapping to df_pattern
df_readR_mapped <- df_pattern %>%
  left_join(df_readR_processed, by = c('seq' = 'seq_target')) %>% 
  dplyr::filter(!is.na(sample_number)) %>% print()



# calculate size factor fo RPM
df_SizeFactor <- df_readR_mapped %>% 
  dplyr::count(sample_number) %>% 
  dplyr::transmute(
    sample_number,
    # calculate size factors for RPM
    sf=1000000/n
  ) %>% 
  print()



# calculate RPM
df_tmp_RPM <- df_readR_mapped %>% 
  # count reads per pattern
  dplyr::group_by(sample_number, pattern_number) %>% 
  dplyr::reframe(reads = n()) %>% 
  dplyr::left_join(df_SizeFactor, by = join_by(sample_number)) %>% 
  dplyr::mutate(RPM = reads * sf) %>% 
  # change RPM to 0 if t1/t2 reads are too low (using cutoff value)
  dplyr::mutate(
    RPM = if_else(
      reads >= read_cutoff, 
      RPM,
      NA
    )
  ) %>% 
  print()


# formatting
df_RPM <- df_pattern %>%
  left_join(
    df_sample_sheet,
    by = 'object',
    relationship = "many-to-many"
  ) %>% 
  dplyr::left_join(df_tmp_RPM, by = join_by(sample_number, pattern_number)) %>% 
  dplyr::select(sample_number, sample_name, org, bc7, bc5, everything()) %>%
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
  ) %>% print()


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
    FC_RPM = log2(t2_cfu/t1_cfu * t2_RPM /t1_RPM ) / (t2_culture_min/60),
    # standardization, Min-Max normalization
    FC_RPM_std = scale(FC_RPM, center = min(FC_RPM, na.rm=T), scale = max(FC_RPM, na.rm=T) - min(FC_RPM, na.rm=T))[,1]
  ) %>% 
  ungroup()


df_change %>% 
  dplyr::select(object, org, pattern_number, wtmt, mt_res, codon, seq, seq_aa, HasStopCodon, antibiotic, conc_antibiotic, rep_cell_lib, rep_selection, rep_DNA_lib,  t1_reads, t2_reads, t1_cfu, t2_cfu, t2_culture_min, t1_RPM, t2_RPM, FC_RPM, FC_RPM_std) %>% 
  dplyr::arrange(pattern_number, rep_selection, conc_antibiotic) %>% 
  write_csv(output_fc)



# optional information ----------------------------------------------------

# output removed raws
df_change %>% 
  filter(is.na(FC_RPM)) %>% 
  dplyr::select(object, org, pattern_number, wtmt, mt_res, codon, seq, seq_aa, HasStopCodon, antibiotic, conc_antibiotic, rep_cell_lib, rep_selection, rep_DNA_lib, t1_cfu, t2_cfu, t2_culture_min, t1_reads, t2_reads, t1_RPM, t2_RPM, FC_RPM) %>% 
  write_csv(output_filter_out)


# output lacked codon/aa information
df_codon_pattern <- df_pattern %>% 
  group_by(seq_aa) %>% 
  dplyr::count(name = 'n_codon_pattern')
  
df_change %>% 
  dplyr::filter(is.na(FC_RPM)) %>% 
  dplyr::select(sample_number, sample_name, seq_aa, codon) %>% 
  dplyr::distinct() %>% 
  add_count(sample_number, sample_name, seq_aa, name = 'n_lacked_codon') %>% 
  left_join(df_codon_pattern) %>% 
  select(-codon) %>% 
  dplyr::distinct() %>% 
  mutate(
    n_in_library = n_codon_pattern - n_lacked_codon,
    InLibrary = case_when(
      n_in_library > 0 ~ 'codon_lacked',
      n_in_library == 0 ~ 'aa_lacked',
    )
    ) %>% 
  tidyr::complete(sample_name, seq_aa) %>% 
  select(sample_number, everything()) %>% 
  arrange(sample_number, seq_aa) %>% 
  write_csv(output_lacked)

# in a subset of the libraries
# in all of the libraries

# normalization (正規化) = Min-Max normalization 
# # FC_RPM_nrm = (FC_RPMa - min(FC_RPMa)) / (max(FC_RPMa) - min(FC_RPMa)),
# # standardization (標準化) = Z-score normalization
# # FC_RPM_std1 = (FC_RPMa - mean(FC_RPMa)) / sd(FC_RPMa),
# # FC_RPM_std2 = scale(FC_RPMa),
# # standardization,
# # FC_RPM_scl1 = scale(FC_RPMa, center = min(FC_RPMa, na.rm=T), scale = max(FC_RPMa, na.rm=T) - min(FC_RPMa, na.rm=T)),
# # FC_RPM_scl2 = scale(FC_RPMa, center = T),
# # # 
