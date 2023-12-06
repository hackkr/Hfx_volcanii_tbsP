library(tidyverse)
library(vroom)

# load in sample key file (also used in mosaics peak calling)
sample_file <- read_csv("03_plot_peaks/03b_hvo_sample_key_extended.csv")

# normalize and calc log2
sample_file %>%
  filter(strain == "WT" & condition == "no.glu") -> WT0
sample_file %>%
  filter(strain == "WT" & condition == "glu") -> WT0.1
sample_file %>%
  filter(strain == "tbsP_HA" & condition == "glu") -> t0.1
sample_file %>%
  filter(strain == "tbsP_HA" & condition == "no.glu") -> t0

## WT no glu
WT.noglu <- left_join(vroom(WT0[1, 4], col_names = F), vroom(WT0[1, 5], col_names = F),
  by = c("X1", "X2")
) %>%
  left_join(vroom(WT0[2, 4], col_names = F), by = c("X1", "X2")) %>%
  left_join(vroom(WT0[2, 5], col_names = F), by = c("X1", "X2")) %>%
  mutate(
    X3.x = .[[3]] + 1, X3.y = .[[4]] + 1, # add 1 read to avoid inf values
    X3.x.x = .[[5]] + 1, X3.y.y = .[[6]] + 1
  ) %>%
  mutate(
    X3.x = .[[3]] / sum(.[[3]]) * 1000000, # normalize each sample by read depth
    X3.y = .[[4]] / sum(.[[4]]) * 1000000,
    X3.x.x = .[[5]] / sum(.[[5]]) * 1000000,
    X3.y.y = .[[6]] / sum(.[[6]]) * 1000000
  ) %>%
  mutate(
    ip_input_A = X3.x / X3.y, # take ratio of IP/WCE
    ip_input_B = X3.x.x / X3.y.y
  ) %>%
  mutate(
    ip_input_mean = rowMeans(.[, 7:8]), # average across bioreps
    counts_log2 = log2(ip_input_mean)
  ) %>% # log2 ratio
  dplyr::select(X1, X2, ip_input_mean, counts_log2) # chrom, position, log2 value

## WT glu
WT.glu <- left_join(vroom(WT0.1[1, 4], col_names = F), vroom(WT0.1[1, 5], col_names = F),
  by = c("X1", "X2")
) %>%
  left_join(vroom(WT0.1[2, 4], col_names = F), by = c("X1", "X2")) %>%
  left_join(vroom(WT0.1[2, 5], col_names = F), by = c("X1", "X2")) %>%
  mutate(
    X3.x = .[[3]] + 1, X3.y = .[[4]] + 1, # add 1 read to avoid inf values
    X3.x.x = .[[5]] + 1, X3.y.y = .[[6]] + 1
  ) %>%
  mutate(
    X3.x = .[[3]] / sum(.[[3]]) * 1000000, # normalize each sample by read depth
    X3.y = .[[4]] / sum(.[[4]]) * 1000000,
    X3.x.x = .[[5]] / sum(.[[5]]) * 1000000,
    X3.y.y = .[[6]] / sum(.[[6]]) * 1000000
  ) %>%
  mutate(
    ip_input_A = X3.x / X3.y, # take ratio of IP/WCE
    ip_input_B = X3.x.x / X3.y.y
  ) %>%
  mutate(
    ip_input_mean = rowMeans(.[, 7:8]), # average across bioreps
    counts_log2 = log2(ip_input_mean)
  ) %>% # log2 ratio
  dplyr::select(X1, X2, counts_log2) # chrom, position, log2 value

## TbsP no glu
t.noglu <- left_join(vroom(t0[1, 4], col_names = F), vroom(t0[1, 5], col_names = F),
  by = c("X1", "X2")
) %>%
  left_join(vroom(t0[2, 4], col_names = F), by = c("X1", "X2")) %>%
  left_join(vroom(t0[2, 5], col_names = F), by = c("X1", "X2")) %>%
  left_join(vroom(t0[3, 4], col_names = F), by = c("X1", "X2")) %>%
  left_join(vroom(t0[3, 5], col_names = F), by = c("X1", "X2")) %>%
  left_join(vroom(t0[4, 4], col_names = F), by = c("X1", "X2")) %>%
  left_join(vroom(t0[4, 5], col_names = F), by = c("X1", "X2")) %>%
  mutate(
    X3.x = .[[3]] + 1, X3.y = .[[4]] + 1, # add 1 read to avoid inf values
    X3.x.x = .[[5]] + 1, X3.y.y = .[[6]] + 1,
    X3.x.x.x = .[[7]] + 1, X3.y.y.y = .[[8]] + 1,
    X3.x.x.x.x = .[[9]] + 1, X3.y.y.y.y = .[[10]] + 1
  ) %>%
  mutate(
    X3.x = .[[3]] / sum(.[[3]]) * 1000000, # normalize each sample by read depth
    X3.y = .[[4]] / sum(.[[4]]) * 1000000,
    X3.x.x = .[[5]] / sum(.[[5]]) * 1000000,
    X3.y.y = .[[6]] / sum(.[[6]]) * 1000000,
    X3.x.x.x = .[[7]] / sum(.[[7]]) * 1000000,
    X3.y.y.y = .[[8]] / sum(.[[8]]) * 1000000,
    X3.x.x.x.x = .[[9]] / sum(.[[9]]) * 1000000,
    X3.y.y.y.y = .[[10]] / sum(.[[10]]) * 1000000
  ) %>%
  mutate(
    ip_input_A = X3.x / X3.y, # take ratio of IP/WCE
    ip_input_B = X3.x.x / X3.y.y,
    ip_input_C = X3.x.x.x / X3.y.y.y,
    ip_input_D = X3.x.x.x.x / X3.y.y.y.y
  ) %>%
  mutate(
    ip_input_mean = rowMeans(.[, 11:14]), # average across bioreps
    counts_log2 = log2(ip_input_mean)
  ) %>% # log2 ratio
  dplyr::select(X1, X2, counts_log2) # chrom, position, log2 value

## TbsP glu
t.glu <- left_join(vroom(t0.1[1, 4], col_names = F), vroom(t0.1[1, 5], col_names = F),
  by = c("X1", "X2")
) %>%
  left_join(vroom(t0.1[2, 4], col_names = F), by = c("X1", "X2")) %>%
  left_join(vroom(t0.1[2, 5], col_names = F), by = c("X1", "X2")) %>%
  left_join(vroom(t0.1[3, 4], col_names = F), by = c("X1", "X2")) %>%
  left_join(vroom(t0.1[3, 5], col_names = F), by = c("X1", "X2")) %>%
  left_join(vroom(t0.1[4, 4], col_names = F), by = c("X1", "X2")) %>%
  left_join(vroom(t0.1[4, 5], col_names = F), by = c("X1", "X2")) %>%
  mutate(
    X3.x = .[[3]] + 1, X3.y = .[[4]] + 1, # add 1 read to avoid inf values
    X3.x.x = .[[5]] + 1, X3.y.y = .[[6]] + 1,
    X3.x.x.x = .[[7]] + 1, X3.y.y.y = .[[8]] + 1,
    X3.x.x.x.x = .[[9]] + 1, X3.y.y.y.y = .[[10]] + 1
  ) %>%
  mutate(
    X3.x = .[[3]] / sum(.[[3]]) * 1000000, # normalize each sample by read depth
    X3.y = .[[4]] / sum(.[[4]]) * 1000000,
    X3.x.x = .[[5]] / sum(.[[5]]) * 1000000,
    X3.y.y = .[[6]] / sum(.[[6]]) * 1000000,
    X3.x.x.x = .[[7]] / sum(.[[7]]) * 1000000,
    X3.y.y.y = .[[8]] / sum(.[[8]]) * 1000000,
    X3.x.x.x.x = .[[9]] / sum(.[[9]]) * 1000000,
    X3.y.y.y.y = .[[10]] / sum(.[[10]]) * 1000000
  ) %>%
  mutate(
    ip_input_A = X3.x / X3.y, # take ratio of IP/WCE
    ip_input_B = X3.x.x / X3.y.y,
    ip_input_C = X3.x.x.x / X3.y.y.y,
    ip_input_D = X3.x.x.x.x / X3.y.y.y.y
  ) %>%
  mutate(
    ip_input_mean = rowMeans(.[, 11:14]), # average across bioreps
    counts_log2 = log2(ip_input_mean)
  ) %>% # log2 ratio
  dplyr::select(X1, X2, counts_log2) # chrom, position, log2 value

# merge into 1 df and export
coverage.all <- data.frame(
  "chr" = WT.noglu$X1,
  "position" = WT.noglu$X2,
  "WT.noglu.counts.log2" = WT.noglu$counts_log2,
  "WT.glu.counts.log2" = WT.glu$counts_log2,
  "tbsP.noglu.counts.log2" = t.noglu$counts_log2,
  "tbsP.glu.counts.log2" = t.glu$counts_log2
)

write_csv(coverage.all, "03_plot_peaks/03b_log2_norm_bedgraph_extended.csv")
