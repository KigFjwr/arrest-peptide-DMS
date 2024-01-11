
# input -------------------------------------------------------------------

args <- commandArgs(trailingOnly = T)
# argument 1 & 2 = read data
input_read1 <- args[1]
input_read2 <- args[2]
# argument 3 = md5 list from source
input_md5 <- args[3]

# output name
output <- 'output/check_md5/check_md5.csv'

# load package ------------------------------------------------------------

library(tidyverse)
library(tools)

# checksum_function -------------------------------------------------------

# prepare output directory
if(!dir.exists('output/')){
  dir.create('output/')
}
if(!dir.exists('output/check_md5/')){
  dir.create('output/check_md5/')
}

df_md5 <- tibble(path = c(input_read1, input_read2)) %>% 
  # get file path info
  dplyr::mutate(
    file = basename(path),
    dir = stringr::str_remove(path, file)
  ) %>% 
  dplyr::select(file, path) %>% 
  
  # calc md5sum
  dplyr::mutate(md5sum_dl = as.character(tools::md5sum(path))) %>% 
  dplyr::left_join(
    readr::read_delim(
      input_md5, 
      col_names = c('md5_source', 'file')
    )
  ) %>% 
  rowwise() %>% 
  dplyr::mutate(
    identical_checksum = identical(md5_source, md5sum_dl),
  ) %>% print()


df_md5 %>% 
  readr::write_csv(output)

