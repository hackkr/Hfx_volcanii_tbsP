TbsP shared peaks (bedgraph)
================
Rylee Hackley

``` r
library(tidyverse)
```

    ## Warning: package 'ggplot2' was built under R version 4.3.1

    ## Warning: package 'purrr' was built under R version 4.3.1

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.2     ✔ readr     2.1.4
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.0
    ## ✔ ggplot2   3.4.3     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.2     ✔ tidyr     1.3.0
    ## ✔ purrr     1.0.1     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(viridis)
```

    ## Warning: package 'viridis' was built under R version 4.3.1

    ## Loading required package: viridisLite

``` r
logs2 <- read_csv("03_plot_peaks/03b_log2_norm_bedgraph_extended.csv")
```

    ## Rows: 4012900 Columns: 6
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (1): chr
    ## dbl (5): position, WT.noglu.counts.log2, WT.glu.counts.log2, tbsP.noglu.coun...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# read in gene boundaries/positions
gff_cds <- read_csv("00_genome_files/genomic_gff_key.csv")
```

    ## Rows: 4066 Columns: 10
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (7): chr, acc, locus_tag, old_locus_tag, strand, type, annotation
    ## dbl (3): length_nt, start, end
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
tbsp <- gff_cds[gff_cds$locus_tag == "HVO_RS02855", ]
gap2 <- gff_cds[gff_cds$locus_tag == "HVO_RS07010", ]
pyrE2 <- gff_cds[gff_cds$locus_tag == "HVO_RS06305", ]
fbp <- gff_cds[gff_cds$locus_tag == "HVO_RS11680", ]
ppsa <- gff_cds[gff_cds$locus_tag == "HVO_RS08595", ]
```

blacklist pyrE2 and CRISPR regions - there’s some weird spikiness around
deletion boundaries and in the CRISPR arrays:

``` r
# pyre2
logs2 %>%
  filter(chr == "NC_013967.1") %>%
  filter(position >= pyrE2$start + 70 & position <= pyrE2$end - 150) %>%
  ggplot() +
  geom_line(aes(
    x = position,
    y = -((tbsP.glu.counts.log2) - (WT.glu.counts.log2)),
    color = "tbsP glu"
  ), size = 0.5, alpha = 0.9) +
  geom_line(aes(
    x = position,
    y = -((tbsP.noglu.counts.log2) - (WT.noglu.counts.log2)),
    color = "tbsP no glu"
  ), size = 0.5) +
  labs(x = "", y = "log2(IP/input)") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_color_manual(
    name = NULL,
    values = c("tbsP no glu" = "#4c1d4b", "tbsP glu" = "#F69c73")
  ) +
  theme_bw()
```

    ## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

![](03c_plot_peaks_bedgraph_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
# CRISPR arrays
logs2 %>%
  filter(chr == "NC_013967.1") %>%
  filter(position >= 2385350 & position <= 2385720) %>%
  ggplot() +
  geom_line(aes(
    x = position,
    y = -((tbsP.glu.counts.log2) - (WT.glu.counts.log2)),
    color = "tbsP glu"
  ), size = 0.5, alpha = 0.9) +
  geom_line(aes(
    x = position,
    y = -((tbsP.noglu.counts.log2) - (WT.noglu.counts.log2)),
    color = "tbsP no glu"
  ), size = 0.5) +
  labs(x = "", y = "log2(IP/input)") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_color_manual(
    name = NULL,
    values = c("tbsP no glu" = "#4c1d4bff", "tbsP glu" = "#F69c73")
  ) +
  theme_bw()
```

![](03c_plot_peaks_bedgraph_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->

``` r
logs2 %>%
  filter(chr == "NC_013967.1") %>%
  filter(position >= 2386360 & position <= 2386500) %>%
  ggplot() +
  geom_line(aes(
    x = position,
    y = -((tbsP.glu.counts.log2) - (WT.glu.counts.log2)),
    color = "tbsP glu"
  ), size = 0.5, alpha = 0.9) +
  geom_line(aes(
    x = position,
    y = -((tbsP.noglu.counts.log2) - (WT.noglu.counts.log2)),
    color = "tbsP no glu"
  ), size = 0.5) +
  labs(x = "", y = "log2(IP/input)") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_color_manual(
    name = NULL,
    values = c("tbsP no glu" = "#4c1d4bff", "tbsP glu" = "#F69c73")
  ) +
  theme_bw()
```

![](03c_plot_peaks_bedgraph_files/figure-gfm/unnamed-chunk-2-3.png)<!-- -->

``` r
logs2 %>%
  filter(chr == "NC_013966.1") %>%
  filter(position >= 204975 + 400 & position <= 207584 - 880) %>%
  ggplot() +
  geom_line(aes(
    x = position,
    y = -((tbsP.glu.counts.log2) - (WT.glu.counts.log2)),
    color = "tbsP glu"
  ), size = 0.5, alpha = 0.9) +
  geom_line(aes(
    x = position,
    y = -((tbsP.noglu.counts.log2) - (WT.noglu.counts.log2)),
    color = "tbsP no glu"
  ), size = 0.5) +
  labs(x = "", y = "log2(IP/input)") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_color_manual(
    name = NULL,
    values = c("tbsP no glu" = "#4c1d4bff", "tbsP glu" = "#F69c73")
  ) +
  theme_bw()
```

![](03c_plot_peaks_bedgraph_files/figure-gfm/unnamed-chunk-2-4.png)<!-- -->

``` r
logs2 %>%
  filter(chr == "NC_013967.1") %>%
  filter(position <= 301800 | position >= 302200) %>%
  filter(position <= 2385250 | position >= 2385720) %>%
  filter(position <= 2386320 | position >= 2386500) -> nc.67

logs2 %>%
  filter(chr == "NC_013966.1") %>%
  filter(position <= 205375 | position >= 206704) -> nc.66

logs.blacklist <- logs2 %>%
  filter(chr != "NC_013966.1" & chr != "NC_013967.1") %>%
  bind_rows(., nc.67, nc.66)

# compare blacklist versus not - 2 peaks removed
logs2 %>%
  filter(chr == "NC_013967.1") %>%
  filter(row_number() %% 20 == 1) %>%
  ggplot() +
  geom_hline(color = "black", yintercept = 2, linetype = "dashed") +
  ylim(c(-1, 2.5)) +
  geom_line(aes(
    x = position,
    y = -(tbsP.glu.counts.log2 - WT.glu.counts.log2),
    color = "tbsP glu"
  ), size = 1.1, alpha = 0.7) +
  geom_line(aes(
    x = position,
    y = -(tbsP.noglu.counts.log2 - WT.noglu.counts.log2),
    color = "tbsP no glu"
  ), size = 0.9, alpha = 0.9) +
  labs(x = "", y = "normalized log2(IP/input)") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_color_manual(
    name = NULL,
    values = c("tbsP no glu" = "#4c1d4bff", "tbsP glu" = "#F69c73")
  ) +
  theme_bw() +
  theme(
    legend.position = c(0.9, .1),
    legend.background = element_rect(color = "white"),
    legend.margin = margin(0, 1, 1, 1)
  )
```

![](03c_plot_peaks_bedgraph_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
logs.blacklist %>%
  filter(chr == "NC_013967.1") %>%
  filter(row_number() %% 20 == 1) %>%
  ggplot() +
  geom_hline(color = "black", yintercept = 2, linetype = "dashed") +
  ylim(c(-1, 2.5)) +
  geom_line(aes(
    x = position,
    y = -(tbsP.glu.counts.log2 - WT.glu.counts.log2),
    color = "tbsP glu"
  ), size = 1.1, alpha = 0.7) +
  geom_line(aes(
    x = position,
    y = -(tbsP.noglu.counts.log2 - WT.noglu.counts.log2),
    color = "tbsP no glu"
  ), size = 0.9, alpha = 0.9) +
  labs(x = "", y = "normalized log2(IP/input)") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_color_manual(
    name = NULL,
    values = c("tbsP no glu" = "#4c1d4bff", "tbsP glu" = "#F69c73")
  ) +
  theme_bw() +
  theme(
    legend.position = c(0.9, .1),
    legend.background = element_rect(color = "white"),
    legend.margin = margin(0, 1, 1, 1)
  )
```

![](03c_plot_peaks_bedgraph_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->

visualize entire chromosomes

``` r
png("03_plot_peaks/NC_013967.1.png", width = 12, height = 4, units = "in", res = 1200)
logs.blacklist %>%
  filter(chr == "NC_013967.1") %>%
  filter(row_number() %% 20 == 1) %>%
  ggplot() +
  geom_hline(color = "black", yintercept = 2, linetype = "dashed") +
  ylim(c(-1, 2.5)) +
  geom_line(aes(
    x = position,
    y = -(tbsP.glu.counts.log2 - WT.glu.counts.log2),
    color = "tbsP glu"
  ), size = 1.1, alpha = 0.7) +
  geom_line(aes(
    x = position,
    y = -(tbsP.noglu.counts.log2 - WT.noglu.counts.log2),
    color = "tbsP no glu"
  ), size = 0.9, alpha = 0.9) +
  labs(x = "", y = "normalized log2(TbsP/WT)") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_color_manual(
    name = NULL,
    values = c("tbsP no glu" = "#4c1d4bff", "tbsP glu" = "#F69c73")
  ) +
  theme_bw() +
  theme(
    legend.position = c(0.9, .1), legend.background = element_rect(color = "white"), legend.margin = margin(0, 1, 1, 1)
  )
dev.off()
```

    ## png 
    ##   2

``` r
png("03_plot_peaks/NC_013966.1.png", width = 8, height = 4, units = "in", res = 1200)
logs.blacklist %>%
  filter(chr == "NC_013966.1") %>%
  filter(row_number() %% 20 == 1) %>%
  ggplot() +
  geom_hline(color = "black", yintercept = 2, linetype = "dashed") +
  ylim(c(-1, 2.5)) +
  geom_line(aes(
    x = position,
    y = -(tbsP.glu.counts.log2 - WT.glu.counts.log2),
    color = "tbsP glu"
  ), size = 1.1, alpha = 0.7) +
  geom_line(aes(
    x = position,
    y = -(tbsP.noglu.counts.log2 - WT.noglu.counts.log2),
    color = "tbsP no glu"
  ), size = 0.9, alpha = 0.9) +
  labs(x = "", y = "normalized log2(TbsP/WT)") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_color_manual(
    name = NULL,
    values = c("tbsP no glu" = "#4c1d4bff", "tbsP glu" = "#F69c73")
  ) +
  theme_bw() +
  theme(
    legend.position = "none"
  )
dev.off()
```

    ## png 
    ##   2

``` r
png("03_plot_peaks/NC_013964.1.png", width = 4, height = 4, units = "in", res = 1200)
logs.blacklist %>%
  filter(chr == "NC_013964.1") %>%
  filter(row_number() %% 20 == 1) %>%
  ggplot() +
  geom_hline(color = "black", yintercept = 2, linetype = "dashed") +
  ylim(c(-1, 2.5)) +
  geom_line(aes(
    x = position,
    y = -(tbsP.glu.counts.log2 - WT.glu.counts.log2),
    color = "tbsP glu"
  ), size = 1.1, alpha = 0.7) +
  geom_line(aes(
    x = position,
    y = -(tbsP.noglu.counts.log2 - WT.noglu.counts.log2),
    color = "tbsP no glu"
  ), size = 0.9, alpha = 0.9) +
  labs(x = "", y = "normalized log2(TbsP/WT)") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_color_manual(
    name = NULL,
    values = c("tbsP no glu" = "#4c1d4bff", "tbsP glu" = "#F69c73")
  ) +
  theme_bw() +
  theme(
    legend.position = "none"
  )
dev.off()
```

    ## png 
    ##   2

do the same for each of the 9 main peaks

``` r
asrna <- rtracklayer::import.bed("00_genome_files/gelsinger2018.bed")
peaks <- read_csv("03_plot_peaks/ChIP_consensus.csv")
```

    ## Rows: 9 Columns: 6
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (2): Chr, strand
    ## dbl (4): start, end, rank, score
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
for (i in 1:nrow(peaks)) {
  filename <- paste("03_plot_peaks/extended_bedgraph_log2/peak_", peaks$rank[i], ".pdf", sep = "")

  gff_cds %>%
    filter(chr == peaks$Chr[i]) %>%
    filter(end >= peaks$start[i] - 1000) %>%
    filter(start <= peaks$end[i] + 1000) -> genes

  asrna %>%
    as.data.frame() %>%
    filter(seqnames == peaks$Chr[i]) %>%
    filter(end >= peaks$start[i] - 1000) %>%
    filter(start <= peaks$end[i] + 1000) -> as.rna

  logs.blacklist %>%
    filter(chr == peaks$Chr[i]) %>%
    filter(position > genes$start[1] & position < tail(genes, n = 1)$end) %>%
    ggplot() +
    geom_line(aes(
      x = position,
      y = -(tbsP.noglu.counts.log2 - WT.noglu.counts.log2),
      color = "tbsP no glu"
    ), size = 1.1) +
    geom_line(aes(
      x = position,
      y = -(tbsP.glu.counts.log2 - WT.glu.counts.log2),
      color = "tbsP glu"
    ), size = 1) +
    labs(x = "", y = "normalized log2(TbsP/WT)") +
    scale_x_continuous(expand = c(0, 0)) +
    ylim(c(-1, 2.5)) +
    geom_hline(yintercept = -0.8, color = "black") +
    scale_color_manual(
      name = NULL,
      values = c("tbsP no glu" = "#4c1d4bff", "tbsP glu" = "#F69c73")
    ) +
    theme_bw() +
    theme(legend.position = c(0.8, 0.75)) -> plot

  # add genes
  plot + geom_rect(
    data = genes,
    aes(xmin = start, xmax = end, ymin = -0.7, ymax = -0.9, fill = strand),
    color = "black", linewidth = 0.5
  ) +
    geom_segment(
      data = as.rna,
      aes(x = start, xend = end, y = -1, yend = -1, color = strand),
      size = 2
    ) -> plot

  pdf(filename, width = 4, height = 4)
  print(plot)
  dev.off()
}
```
