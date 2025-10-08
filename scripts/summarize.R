# last update: 2025/10/01

# input -------------------------------------------------------------------

args <- commandArgs(trailingOnly = T)

default_args <- c(
  "", 
  "", 
  "", 
  "", 
  "", 
  "", 
  "conc_antibiotic",
  "0",
  'output/fig/reproducibility_growth_rate.pdf',
  'output/fig/reproducibility_fitness.pdf',
  'output/fig/heatmap_growthrate.pdf',
  'output/fig/heatmap_fitness.pdf',
  'input/plot_order.csv'
)   # default setting
default_flg <- is.na(args[1:13])
args[default_flg] <- default_args[default_flg]  

# argument 1: target name
# example: input_target <- 'ApdA_S109toP123_PY79'
input_target <- args[1]

# argument 2: DMS_start residue number
# example: input_target_AA_number_s <- 109
input_target_AA_number_s <- as.numeric(args[2])

# argument 3: DMS_end residue number
# example: input_target_AA_number_e <- 123
input_target_AA_number_e <- as.numeric(args[3])

# argument 4: input csv file, fitness list
# example: input_calc_GR <- 'output/calc/CalcGR_ApdA_S109-P123_PY79_qc30_co10.csv'
input_calc_GR <- args[4]

# argument 5: input csv file, mutation pattern list
# example: input_pattern <- 'output/ptn_ApdA_S109toP123_NNK.csv'
input_pattern <-  args[5]

# argument 6: output suffix
# example: output_suffix <- 'qc30_co10'
output_suffix <-  args[6]


# optional arguments

# argument 7: baseline_column
# example: baseline_col <- "conc_antibiotic"
baseline_col <- args[7]

# argument 8: baseline_value
# example: baseline_value <- 0
baseline_value <- as.numeric(args[8])

# argument 9: output plot file name, scatter plot
# example: output_plot_rep_gr <- 'output/fig/reproducibility_growth_rate.pdf'
output_plot_rep_gr <- args[9]

# argument 10: output plot file name, scatter plot
# example: output_plot_rep_fn <- 'output/fig/reproducibility_fitness.pdf'
output_plot_rep_fn <- args[10]

# argument 11: output plot file name, heat map
# example: output_plot_hm_gr <- 'output/fig/heatmap_growthrate.pdf'
output_plot_hm_gr <- args[11]

# argument 12: output plot file name, heat map
# example: output_plot_hm_fn <- 'output/fig/heatmap_fitness.pdf'
output_plot_hm_fn <- args[12]

# argument 13: input file name, plot order
# example: input_plot_order <- 'input/plot_order.csv'
input_plot_order <- args[13]


# load package ------------------------------------------------------------

library(tidyverse)
library(Biostrings)
library(ggthemes)
library(patchwork)

# source ------------------------------------------------------------------



# function ----------------------------------------------------------------



# scatter plot to check reproducibility

fx_plot_reproducibility <- function(x_df, y_title, z_subtitle){
  xy_min <- min(c(x_df$rep_1, x_df$rep_2), na.rm = TRUE)
  xy_max <- max(c(x_df$rep_1, x_df$rep_2), na.rm = TRUE)
  x_df %>%
    ggplot(aes(rep_1, rep_2)) +
    geom_point(color = 'grey20', alpha = 0.5, shape = 16, size = 0.8) +
    scale_x_continuous(limits = c(xy_min - 0.1, xy_max + 0.1)) +
    scale_y_continuous(limits = c(xy_min - 0.1, xy_max + 0.1)) +
    facet_wrap(~wrap) +
    coord_fixed() +
    ggtitle(y_title, z_subtitle) +
    labs(x='replicate 1', y='replicate 2') +
    theme_few(base_family="Helvetica")
}

# for plotting of heatmap

fx_plot_heatmap <- function(x_df, y_value, za_label, zb_label){
  title1 <- x_df %>% select({{ za_label}}) %>% pull() %>% unique() %>% as.character()
  title2 <- zb_label
  output <- ggplot(x_df, aes(res_plot, reorder( {{y_value}}, -plot_order), fill = mean)) +
    geom_tile() +
    scale_fill_viridis_c()  +
    labs(x = 'position', y = 'replaced to') +
    scale_x_discrete(labels = c(v_label_position, 'wt') %>% unlist) +
    ggtitle(title1, title2) +
    theme_few(base_family="Helvetica") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}






# load input --------------------------------------------------------------



df_GR <- read_csv(input_calc_GR) %>% 
  mutate(
    sel_condition = str_c(antibiotic, '_', conc_antibiotic)
  )





# prepare output directory
if(!dir.exists('output/fig/')){
  dir.create('output/fig/')
}


# prepare output file name
output_plot_rep_gr <- output_plot_rep_gr %>% str_replace('.pdf', str_c('_', output_suffix, '.pdf'))
output_plot_rep_fn <- output_plot_rep_fn %>% str_replace('.pdf', str_c('_', output_suffix, '.pdf'))
output_plot_hm_gr <- output_plot_hm_gr %>% str_replace('.pdf', str_c('_', output_suffix, '.pdf'))
output_plot_hm_fn <- output_plot_hm_fn %>% str_replace('.pdf', str_c('_', output_suffix, '.pdf'))


# formatting --------------------------------------------------------------

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





# prepare amino acid residues with their residue number

input_target_AA <- df_GR %>% filter(wtmt == 'wt') %>% pull(seq_aa) %>% unique()
v_input_target_AA <- strsplit(input_target_AA, '') %>% unlist
v_label_position <- str_c(v_input_target_AA, c(input_target_AA_number_s:input_target_AA_number_e))





# calculate fitness (vs base condition)

# codon level
df_control <- df_GR %>%
  dplyr::select(condition = all_of(baseline_col), everything()) %>% 
  filter(condition == baseline_value) %>% 
  dplyr::select(pattern_number, seq, seq_aa, wtmt, mt_res, codon, rep_selection, value_control = Growth_Rate, condition_control = sel_condition)

df_treated <- df_GR %>% 
  dplyr::select(condition = all_of(baseline_col), everything()) %>% 
  filter(condition != baseline_value) %>% 
  dplyr::select(pattern_number, seq, seq_aa, wtmt, mt_res, codon, rep_selection, value_treated = Growth_Rate, condition_treated = sel_condition)

df_fitness <- df_control %>% 
  right_join(df_treated, by = join_by(pattern_number, seq, seq_aa, wtmt, mt_res, codon, rep_selection)) %>% 
  dplyr::mutate(
    rep_selection=str_c('rep_', rep_selection),
    fitness = value_treated / value_control
  ) %>% 
  dplyr::select(-value_control, -value_treated) %>% print()



# plotting, reproducibility -----------------------------------------------

# Growth rate

# codon level
p_rep_GR_c <- df_GR %>% 
  dplyr::select(pattern_number, seq_aa, wtmt, HasStopCodon, sel_condition, rep_selection, Growth_Rate) %>% 
  dplyr::mutate(rep_selection=str_c('rep_', rep_selection)) %>% 
  tidyr::pivot_wider(names_from = rep_selection, values_from = Growth_Rate) %>% 
  filter(!(is.na(rep_1)|is.na(rep_2))) %>% 
  mutate(wrap = sel_condition) %>% 
  fx_plot_reproducibility('codon level', 'growth rate (GR)')

# aa level
p_rep_GR_a <- df_GR %>% 
  dplyr::select(pattern_number, seq_aa, sel_condition, rep_selection, Growth_Rate) %>% 
  dplyr::mutate(rep_selection=str_c('rep_', rep_selection)) %>% 
  dplyr::group_by(sel_condition, rep_selection, seq_aa) %>% 
  dplyr::reframe(Growth_Rate_aa = mean(Growth_Rate), n=n()) %>% 
  tidyr::pivot_wider(names_from = rep_selection, values_from = Growth_Rate_aa) %>% 
  filter(!(is.na(rep_1)|is.na(rep_2))) %>% 
  mutate(wrap = sel_condition) %>% 
  fx_plot_reproducibility('AA level', 'growth rate (GR)')

p_rep_GR <- p_rep_GR_c/p_rep_GR_a +
  plot_annotation(
    title = input_target,
    subtitle = 'reproducibility',
    caption = Sys.Date()
  )

ggsave(output_plot_rep_gr, p_rep_GR, width = 5, height = 6)



# Fitness
p_rep_FN_c <- df_fitness %>%
  filter(!is.na(fitness)) %>% 
  tidyr::pivot_wider(names_from = rep_selection, values_from = fitness) %>%
  filter(!(is.na(rep_1)|is.na(rep_2))) %>% 
  mutate(wrap = condition_treated) %>%
  fx_plot_reproducibility('codon level', 'GR_treated / GR_control')

p_rep_FN_a <- df_fitness %>%
  filter(!is.na(fitness)) %>% 
  dplyr::group_by(condition_treated, mt_res, rep_selection, seq_aa) %>%
  dplyr::reframe(fitness_aa = mean(fitness, na.rm = T), n=n()) %>%
  tidyr::pivot_wider(names_from = rep_selection, values_from = fitness_aa) %>%
  filter(!(is.na(rep_1)|is.na(rep_2))) %>% 
  mutate(wrap = condition_treated) %>%
  fx_plot_reproducibility('AA level', 'GR_treated / GR_control')

p_rep_FN <- p_rep_FN_c/p_rep_FN_a +
  plot_annotation(
    title = input_target,
    subtitle = 'reproducibility',
    caption = Sys.Date()
  )

ggsave(output_plot_rep_fn, p_rep_FN, width = 5, height = 6)





# plotting, heatmap -------------------------------------------------------



## growth rate
### prepare dataframe for heatmap
### plot
df_plot_heatmap_gr <- df_GR %>% 
    dplyr::select(pattern_number, seq, seq_aa, antibiotic, wtmt, mt_res, codon, sel_condition, rep_selection, Growth_Rate) %>%
    dplyr::mutate(
      rep_selection=str_c('rep_', rep_selection)
    ) %>%
    tidyr::pivot_wider(names_from = rep_selection, values_from = Growth_Rate) %>%
    dplyr::mutate(
      mean = (rep_1+rep_2)/2,
      res = str_extract(mt_res, '\\d+$'),
      # codon = case_when(
      #   mt_res == 'wt' ~ 'wt',
      #   str_sub(mt_res, 1,1) == 'r' ~ str_sub(seq, as.integer(str_extract(mt_res, '\\d+$'))*3-2, as.integer(str_extract(mt_res, '\\d+$'))*3)
      # ),
      res_plot = if_else(
        !is.na(res),
        as.character(input_target_AA_number_s + as.integer(res) - 1),
        res
      )
    ) %>%
    dplyr::left_join(df_ac_order, by = 'codon')

### plot
p_hm_gr_cdn <- df_plot_heatmap_gr %>% fx_plot_heatmap(aa_codon, antibiotic, 'GR, empty: wild-type, grey: missing') + facet_wrap(~sel_condition)
p_hm_gr_aa <- df_plot_heatmap_gr %>% fx_plot_heatmap(aa, antibiotic, 'GR, empty: wild-type, grey: missing') + facet_wrap(~sel_condition)
p_hm_gr <- patchwork::wrap_plots(p_hm_gr_cdn)/patchwork::wrap_plots(p_hm_gr_aa) +
  plot_annotation(
    title = input_target,
    subtitle = 'Heatmap, growth rate (GR)',
    caption = Sys.Date()
  )

### output
ggsave(
  output_plot_hm_gr, 
  p_hm_gr, 
  width = 9*(input_target_AA_number_e - input_target_AA_number_s)/8,
  height = 10
)



## fitness
### prepare dataframe for heatmap
df_plot_heatmap_fn <- df_fitness %>% 
  tidyr::pivot_wider(names_from = rep_selection, values_from = fitness) %>%
  dplyr::mutate(
    mean = (rep_1+rep_2)/2,
    res = str_extract(mt_res, '\\d+$'), 
    res_plot = if_else(
      !is.na(res),
      as.character(input_target_AA_number_s + as.integer(res) - 1),
      res
    )
  ) %>% 
  dplyr::left_join(df_ac_order, by = 'codon') %>%
  group_split(condition_treated)

### plot
p_hm_fn_cdn <- df_plot_heatmap_fn %>% purrr::map(~ fx_plot_heatmap(.x, aa_codon, condition_treated, 'GR_treated / GR_control, empty: wild-type, grey: missing'))
p_hm_fn_aa <- df_plot_heatmap_fn %>% purrr::map(~ fx_plot_heatmap(.x, aa, condition_treated, 'GR_treated / GR_control, empty: wild-type, grey: missing'))
p_hm_fn <- patchwork::wrap_plots(p_hm_fn_cdn)/patchwork::wrap_plots(p_hm_fn_aa)  +
  plot_annotation(
    title = input_target,
    subtitle = 'heatmap, fitness',
    caption = Sys.Date()
  )

### output
ggsave(
  output_plot_hm_fn, 
  p_hm_fn, 
  width = 9*(input_target_AA_number_e - input_target_AA_number_s)/8,
  height = 10
  )

 

# output source data ------------------------------------------------------


df_fitness %>% 
  arrange(rep_selection) %>% 
  write_csv(
    output_plot_hm_fn %>% 
      str_replace('fitness', 'fitness_source') %>% 
      str_replace('pdf', 'csv')
  )  
