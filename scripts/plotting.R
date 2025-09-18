
# input -------------------------------------------------------------------

args <- commandArgs(trailingOnly = T)

default_args <- c("", "", "", "", "", "output/fig/reproducibility.pdf", 'output/fig/heatmap.pdf', 'input/plot_order.csv')   # default setting
default_flg <- is.na(args[1:8])
args[default_flg] <- default_args[default_flg]  

# argument 1: target name, like 'ApdP_Q126-P134_JM109'
input_target <- args[1]

# argument 2: DMS_start residue number, like 126
input_target_AA_number_s <- as.numeric(args[2])

# argument 3: DMS_end residue number, like 134
input_target_AA_number_e <- as.numeric(args[3])

# argument 4: input csv file, fitness list, like 'output/calc/CalcFC_ApdP_Q126-P134_JM109.csv'
input_calc_FC <- args[4]

# argument 5: input csv file, mutation pattern list, like 'output/ptn_ApdP_Q126toP134_NNN.csv'
input_pattern <-  args[5]


# optional arguments
# argument 6: output plot file name, scatter plot
output_plot_rep <- args[6]
# argument 7: output plot file name, heat map
output_plot_heat <- args[7]
# argument 8: input file name, plot order
input_plot_order <- args[8]

# output_plot_het_cd <- 'output/fig/heatmap_codon.pdf'
# output_plot_het_aa <- 'output/fig/heatmap_aa.pdf'

# load package ------------------------------------------------------------

library(tidyverse)
library(Biostrings)
library(patchwork)

# source ------------------------------------------------------------------

# load input --------------------------------------------------------------

df_FC <- read_csv(input_calc_FC) %>% 
  mutate(motif = str_sub(seq_aa, -4,-1))

df_plot <- df_FC %>% 
  filter(!is.na(FC_RPM_std))

df_pattern <- read_csv(input_pattern) %>% 
  mutate(motif = str_sub(seq_aa, -4,-1)) 

# prepare output directory
if(!dir.exists('output/fig/')){
  dir.create('output/fig/')
}

# formatting --------------------------------------------------------------

# df_aa_pattern <- df_pattern %>% 
#   select(object, wtmt, seq_aa, stop, motif) %>%
#   distinct() %>% 
#   filter(!(wtmt=='mt' & seq_aa == 'QSKCIRAPP'))

df_aa_pattern <- df_pattern %>% 
  select(seq_aa) %>% 
  distinct() %>% 
  mutate(aa_pattern_number = row_number()) %>% 
  left_join(
    df_pattern %>% 
      select(seq_aa) %>% 
      distinct() %>% 
      mutate(aa_pattern_number = row_number()) %>% 
      tidyr::separate(
        col = seq_aa,
        into = paste0("r", 0:9),  # 新しい列の名前
        sep = ""
      )
  ) %>% 
  dplyr::select(aa_pattern_number, seq_aa, r1:r9)



# prepare list of amino acid with their order on plot
df_order <- read_csv(input_plot_order)
df_aa_codon <- tibble(
  aa=as.character(Biostrings::GENETIC_CODE),
  codon=as.character(names(Biostrings::GENETIC_CODE))
) %>% 
  dplyr::mutate(aa_codon = str_c(aa, "/",codon)) %>% 
  arrange(aa_codon)


df_aa_order <- df_order %>% select(plot_order, aa=custom_1) %>% filter(!is.na(aa))
df_ac_order <- df_order %>% select(plot_order, aa=custom_1) %>% filter(!is.na(aa)) %>% 
  dplyr::left_join(df_aa_codon) %>% 
  dplyr::mutate(
    aa_codon = if_else(is.na(aa_codon), aa, aa_codon),
    codon = if_else(is.na(codon), aa, codon),
  ) %>% 
  arrange(plot_order, aa, codon) %>% 
  dplyr::mutate(plot_order = row_number()) 


# df_cd_order <- df_order %>% select(plot_order, codon=custom_4) %>% filter(!is.na(codon))




# plotting, reproducibility -----------------------------------------------

# codon
p_rep_c <- df_plot %>% 
  dplyr::select(pattern_number, seq_aa, wtmt, HasStopCodon, motif, conc_antibiotic, rep_selection, FC_RPM_std) %>% 
  dplyr::mutate(rep_selection=str_c('n', rep_selection)) %>% 
  tidyr::pivot_wider(names_from = rep_selection, values_from = FC_RPM_std) %>% 
  dplyr::rename(rep1=n1, rep2=n2) %>% 
  ggplot(aes(rep1, rep2)) +
  geom_point(color = 'grey20', alpha = 0.5, shape = 16, size = 0.8) +
  facet_wrap(~conc_antibiotic) +
  ggtitle('Codon level', 'FC_RPM values (scaled) of each mutant are plotted')


# aa
p_rep_a <- df_plot %>% 
  dplyr::select(pattern_number, seq_aa, conc_antibiotic, rep_selection, FC_RPM_std) %>% 
  dplyr::mutate(rep_selection=str_c('n', rep_selection)) %>% 
  dplyr::group_by(conc_antibiotic, rep_selection, seq_aa) %>% 
  dplyr::reframe(FC_RPM_aa = mean(FC_RPM_std), n=n()) %>% 
  tidyr::pivot_wider(names_from = rep_selection, values_from = FC_RPM_aa) %>% 
  dplyr::rename(rep1=n1, rep2=n2) %>% 
  ggplot(aes(rep1, rep2)) +
  geom_point(color = 'grey20', alpha = 0.5, shape = 16, size = 0.8) +
  facet_wrap(~conc_antibiotic) +
  ggtitle('AA level', 'mean FC_RPM values (scaled) of each amino acid sequence are plotted')

p_rep <- p_rep_c/p_rep_a +
  plot_annotation(
    title = input_target,
    subtitle = 'Reproducibility check',
    caption = Sys.Date()
  )

ggsave(output_plot_rep, p_rep, width = 6, height = 6)


# heatmap -----------------------------------------------------------------

p_heat_codon <- df_FC %>% 
  dplyr::select(seq, mt_res, conc_antibiotic, rep_selection, FC_RPM_std) %>% 
  dplyr::mutate(
    rep_selection=str_c('n', rep_selection)
  ) %>% 
  tidyr::pivot_wider(names_from = rep_selection, values_from = FC_RPM_std) %>% 
  dplyr::mutate(
    mean = (n1+n2)/2,
    res = str_extract(mt_res, '\\d+$'), 
    codon = case_when(
      mt_res == 'wt' ~ 'wt',
      str_sub(mt_res, 1,1) == 'r' ~ str_sub(seq, as.integer(str_extract(mt_res, '\\d+$'))*3-2, as.integer(str_extract(mt_res, '\\d+$'))*3)
    ),
    res_plot = if_else(
      !is.na(res),
      as.character(input_target_AA_number_s + as.integer(res) - 1),
      res
    ),
  ) %>% 
  dplyr::left_join(df_ac_order, by = 'codon') %>%
  ggplot(aes(res_plot, reorder(aa_codon, -plot_order), fill = mean)) +
  geom_tile() +
  scale_fill_viridis_c()  +
  facet_wrap(~conc_antibiotic) +
  labs(x = 'mutated residue number', y = 'replaced to') +
  ggtitle('Codon level', 'mean FC_RPM (scaled) of each codons, empty: wild-type, grey: lacked')



p_heat_aa <- df_FC %>% 
  dplyr::select(seq_aa, mt_res, conc_antibiotic, rep_selection, FC_RPM_std) %>% 
  dplyr::mutate(
    rep_selection=str_c('n', rep_selection),
    aa_pattern_number = row_number()
    ) %>% 
  dplyr::group_by(conc_antibiotic, mt_res, rep_selection, seq_aa) %>% 
  dplyr::reframe(FC_RPM_aa = mean(FC_RPM_std, na.rm = T), n=n()) %>% 
  tidyr::pivot_wider(names_from = rep_selection, values_from = FC_RPM_aa) %>% 
  mutate(
    mean = (n1+n2)/2,
    res = str_extract(mt_res, '\\d+$'), 
    aa = case_when(
      mt_res == 'wt' ~ 'wt',
      str_sub(mt_res, 1,1) == 'r' ~ str_sub(seq_aa, as.integer(str_extract(mt_res, '\\d+$')), as.integer(str_extract(mt_res, '\\d+$')))
    ),
    res_plot = if_else(
      !is.na(res),
      as.character(input_target_AA_number_s + as.integer(res) - 1),
      res
    ),
  ) %>% 
  dplyr::left_join(df_aa_order, by = 'aa') %>% 
  ggplot(aes(res_plot, reorder(aa, -plot_order), fill = mean)) +
  geom_tile() +
  scale_fill_viridis_c() +
  facet_wrap(~conc_antibiotic) +
  labs(x = 'mutated residue number', y = 'replaced to') +
  ggtitle('AA level', 'mean FC_RPM (scaled) of each aa, grey: lacked')

p_heat <- p_heat_codon/p_heat_aa +
  plot_annotation(
    title = input_target,
    subtitle = 'Heatmap',
    caption = Sys.Date()
  ) +
  plot_layout(
    heights = c(2.8,1),
  )

ggsave(
  output_plot_heat, 
  p_heat, 
  width = 9*(input_target_AA_number_e - input_target_AA_number_s)/8, 
  height = 12
)

# ggsave(
#   output_plot_het_cd, 
#   p_heat_codon, 
#   width = 9*(input_target_AA_number_e - input_target_AA_number_s)/10, 
#   height = 9
# )
# ggsave(
#   output_plot_het_aa, 
#   p_heat_aa, 
#   width = 9*(input_target_AA_number_e - input_target_AA_number_s)/10, 
#   height = 6
# )




# Rx ----------------------------------------------------------------------


# 
# df_plot %>% 
#   mutate(
#     color = case_when(
#       wtmt == 'wt' ~ '1',
#       stop == T ~ '2',
#       #wtmt == 'mt' ~ '#FF4B00'
#     )
#   ) %>% 
#   group_by(pattern_number, conc_antibiotic, rep_cell_lib) %>% 
#   mutate(mean_FC = mean(FC_RPM_std)) %>% 
#   ungroup() %>% 
#   #mutate(motif = if_else(stop == T, 'stop', motif)) %>%
#   mutate(motif = if_else(str_detect(seq_aa, 'QSKCIRAPP'), 'WT', motif)) %>%
#   #filter(motif %in% c('RAPP','RGPP', 'AAPP', 'RAPA', 'RAAP', 'RAGP', 'stop', 'WT')) %>%
#   
#   ggplot(aes(x = motif, y = mean_FC, color = color, fill = color)) +
#   geom_point() +
#   # geom_point(
#   #shape = 1,
#   #   size = 1,
#   #   stroke = 1,
#   #   alpha = 0.4,
#   #   position = position_jitter(
#   #     width = 0.1,
#   #     height = 0,
#   #     seed = 328
#   #   )
#   # ) +
# #scale_y_continuous(limits = c(0,2.5))+
# facet_wrap(rep_selection~conc_antibiotic)
# 
# 
# 
# 
# 
# 
# 
# df_plot %>% 
#   select(seq, seq_aa, wtmt, conc_antibiotic, rep_selection, FC_RPM_std) %>% 
#   mutate(rep_selection=str_c('n', rep_selection)) %>% 
#   pivot_wider(names_from = rep_selection, values_from = FC_RPM_std) %>% 
#   filter(conc_antibiotic != 0) %>% 
#   mutate(tmp=n1-n2) %>% 
#   arrange(tmp)
# 
# df_plot %>% 
#   filter(is.na(FC_RPMa)) %>% 
#   select(object:motif) %>% 
#   write_csv('c.csv')
