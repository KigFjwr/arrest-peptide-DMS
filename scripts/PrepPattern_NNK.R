
# input -------------------------------------------------------------------

args <- commandArgs(trailingOnly = T)

# argument 1: input_peptide
# example: input_peptide <- 'CliM'
input_peptide <- args[1]

# argument 2: input sequence which is the target of mutation
# example: input_string <- "AAGTATGTTTTAATAAGAGACATATTTGTAAATAGAATTACATATTCTGAAGAACGACTGCCTAAACAGTATATAGTTTTTCAGAAATATGATATTTGGCGGTATTGTAGTTTATTTAAA"
input_string <- args[2]

# argument 3: output file name
# example: output <- 'output/ptn_CliM_38to77_NNK.csv'
output <- args[3]


# package -----------------------------------------------------------------


library(tidyverse)
library(Biostrings)
library(gtools)


# wrangling ---------------------------------------------------------------


# prepare output directory
if(!dir.exists('output/')){
  dir.create('output/')
}
if(!dir.exists('output/')){
  dir.create('output/')
}




# prep wild type sequence

v_split_codon <- str_sub(input_string, start = seq(1, nchar(input_string), by = 3), end = seq(3, nchar(input_string), by = 3)) %>% print()
len_aa <- length(v_split_codon)

df_wt <- tibble(codon = v_split_codon) %>% 
  dplyr::mutate(
    residues = str_c('r', row_number())
  ) %>% 
  pivot_wider(
    names_from = residues,
    values_from = codon
  ) %>% 
  mutate(
    object = input_peptide,
    wtmt = 'wt',
    codon = 'wt',
    mt_res = 'wt',
    pattern_number = row_number()
  )



# prep codon list
df_all_codon <- expand_grid(
  c1 = c("A", "C", "G", "T"),
  c2 = c("A", "C", "G", "T"),
  c3 = c("G", "T")
) %>%
  transmute(codon = str_c(c1,c2,c3)) %>% 
  dplyr::mutate(wtmt = 'mt') %>%
  dplyr::left_join(
    tibble(
      wtmt = 'mt',
      mt_res = c(1:len_aa)
    ),
    relationship = "many-to-many"
  ) %>% 
  dplyr::arrange(mt_res) %>% 
  rowwise() %>% 
  mutate(
    r = list(replace(v_split_codon, mt_res, codon))
  ) %>% 
  unnest_wider(r, names_sep = '') %>% 
  mutate(mt_res = str_c('r', mt_res))


# 結合するcodonの列名をベクトルで指定する
columns_to_combine <- str_c("r", 1:len_aa)


# add mutant pattern
df_all_pattern <- df_wt %>% 
  
  # mt情報を追加
  dplyr::bind_rows(df_all_codon) %>% 
  
  # order列を調整
  dplyr::mutate(pattern_number = row_number()) %>% 
  
  # object列を埋める
  tidyr::fill(object) %>% 
  
  # columns_to_combineで指定した列の値ををstr_cで結合する
  dplyr::mutate(seq = do.call(stringr::str_c, c(across(all_of(columns_to_combine)), sep = ""))) %>% 
  
  # BiostringsでseqのDNA配列を翻訳し、bind_colsで結合する
  dplyr::bind_cols(
    seq_aa = as.character(
      Biostrings::DNAStringSet(.$seq) %>% 
        Biostrings::translate() 
    )
  ) %>% 
  
  # HasStop列を追加する。stopコドンが現れる場合にはTRUE、ない場合にはFALSE
  dplyr::mutate(HasStopCodon = if_else(str_detect(seq_aa, stringr::fixed('*')), TRUE, FALSE)) %>% 
  
  # WTと同じ配列なのに'mt'となっている行を削除する
  filter(!(seq == input_string & wtmt == 'mt')) %>% 
  print()

df_all_pattern %>% 
  write_csv(output)
  
  
