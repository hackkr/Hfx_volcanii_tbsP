02b_motility
================
Rylee Hackley

``` r
library(tidyverse)

sch <- read_csv("02_schiller_overlap/02_schiller2023_list.csv")

gff.df <- read_csv("000_20230911_hvo/genomic_gff_key.csv")
trmB <- read_csv("01_deseq2_output/trmB_interaction_all.csv") %>%
  filter(padj <= 0.01) %>%
  mutate(model = "TrmB interaction")
tbsP <- read_csv("01_deseq2_output/tbsP_interaction_all.csv") %>%
  filter(padj <= 0.01) %>%
  mutate(model = "TbsP interaction")
tbsP_gen <- read_csv("01_deseq2_output/tbsP_genotype_all.csv") %>%
  filter(padj <= 0.01) %>%
  mutate(model = "TbsP genotype")
tbsP_degs <- full_join(tbsP[c(1, 9, 11)], tbsP_gen[c(1, 9, 11)], by = c("locus_tag", "old_locus_tag"))

sch <- left_join(sch, gff.df[c(3, 4)])
```

``` r
inner_join(tbsP_degs, sch) -> tP
```

    ## Joining with `by = join_by(locus_tag, old_locus_tag)`

``` r
write_csv(tP, "02_schiller_overlap/tbsP_overlap.csv")
inner_join(trmB[c(1, 9, 11)], sch) -> tB
```

    ## Joining with `by = join_by(locus_tag, old_locus_tag)`

``` r
write_csv(tB, "02_schiller_overlap/trmB_overlap.csv")
```

``` r
tP
```

    ## # A tibble: 11 × 9
    ##    locus_tag   old_locus_tag model.x     model.y analysis_type.x analysis_type.y
    ##    <chr>       <chr>         <chr>       <chr>   <chr>           <chr>          
    ##  1 HVO_RS08595 HVO_0812      TbsP inter… <NA>    <NA>            <NA>           
    ##  2 HVO_RS08985 HVO_0894      TbsP inter… <NA>    <NA>            rdfA_ddfA      
    ##  3 HVO_RS19100 HVO_2976      TbsP inter… <NA>    <NA>            <NA>           
    ##  4 HVO_RS01000 HVO_B0203     TbsP inter… <NA>    <NA>            <NA>           
    ##  5 HVO_RS01050 HVO_B0213     TbsP inter… <NA>    <NA>            <NA>           
    ##  6 HVO_RS06840 HVO_0445      <NA>        TbsP g… <NA>            <NA>           
    ##  7 HVO_RS06845 HVO_0446      <NA>        TbsP g… <NA>            <NA>           
    ##  8 HVO_RS10425 HVO_1191      <NA>        TbsP g… <NA>            rdfA_ddfA      
    ##  9 HVO_RS14465 HVO_2031      <NA>        TbsP g… <NA>            rdfA_ddfA      
    ## 10 HVO_RS16065 HVO_2361      <NA>        TbsP g… JK3_cetZ1       <NA>           
    ## 11 HVO_RS18640 HVO_2882      <NA>        TbsP g… <NA>            <NA>           
    ## # ℹ 3 more variables: analysis_type <chr>, uniprot <chr>, annotation <chr>

``` r
tbsP[tbsP$locus_tag %in% tP$locus_tag, ] # genes potentially involved in -glucose hypermotility: ppsA, acs (0894), 2031 implicated in sugar metabolism. AlgJ not included in schiller data. nothing near peaks... what's the effect in schiller? whats the direction of expression change in our data?
```

    ## # A tibble: 5 × 11
    ##   locus_tag   baseMean log2FoldChange lfcSE  stat       pvalue      padj acc    
    ##   <chr>          <dbl>          <dbl> <dbl> <dbl>        <dbl>     <dbl> <chr>  
    ## 1 HVO_RS08595   40407.           2.22 0.241  5.07 0.000000395  0.000105  WP_004…
    ## 2 HVO_RS08985   12133.          -3.01 0.471 -4.26 0.0000201    0.00291   WP_004…
    ## 3 HVO_RS19100    4296.          -3.07 0.490 -4.23 0.0000233    0.00316   WP_004…
    ## 4 HVO_RS01000    1346.          -3.39 0.444 -5.39 0.0000000714 0.0000272 WP_013…
    ## 5 HVO_RS01050    6475.          -3.61 0.537 -4.85 0.00000121   0.000252  WP_013…
    ## # ℹ 3 more variables: old_locus_tag <chr>, annotation <chr>, model <chr>

``` r
tbsP_gen[tbsP_gen$locus_tag %in% tP$locus_tag, ]
```

    ## # A tibble: 6 × 11
    ##   locus_tag   baseMean log2FoldChange lfcSE   stat   pvalue     padj acc        
    ##   <chr>          <dbl>          <dbl> <dbl>  <dbl>    <dbl>    <dbl> <chr>      
    ## 1 HVO_RS06840   10894.           2.43 0.297   4.82 1.43e- 6 1.96e- 4 WP_0040444…
    ## 2 HVO_RS06845   11249.           2.50 0.315   4.76 1.91e- 6 2.54e- 4 WP_0040444…
    ## 3 HVO_RS10425   13621.          -5.16 0.380 -11.0  6.62e-28 4.39e-25 WP_0040437…
    ## 4 HVO_RS14465    2094.          -3.83 0.452  -6.27 3.62e-10 6.55e- 8 WP_0040419…
    ## 5 HVO_RS16065   35078.           2.07 0.177   6.04 1.53e- 9 2.54e- 7 WP_0040422…
    ## 6 HVO_RS18640    8597.           2.12 0.115   9.75 1.93e-22 6.40e-20 WP_0130355…
    ## # ℹ 3 more variables: old_locus_tag <chr>, annotation <chr>, model <chr>

``` r
tB
```

    ## # A tibble: 28 × 8
    ##    locus_tag   old_locus_tag model analysis_type.x analysis_type.y analysis_type
    ##    <chr>       <chr>         <chr> <chr>           <chr>           <chr>        
    ##  1 HVO_RS05875 HVO_0240      TrmB… <NA>            <NA>            clusters_14  
    ##  2 HVO_RS05875 HVO_0240      TrmB… <NA>            <NA>            clusters_12  
    ##  3 HVO_RS06270 HVO_0326      TrmB… <NA>            <NA>            clusters_12  
    ##  4 HVO_RS06310 HVO_0334      TrmB… <NA>            rdfA_ddfA       <NA>         
    ##  5 HVO_RS06430 HVO_0358      TrmB… <NA>            <NA>            clusters_14  
    ##  6 HVO_RS08595 HVO_0812      TrmB… <NA>            <NA>            clusters_14  
    ##  7 HVO_RS09725 HVO_1048      TrmB… JK3_cetZ1       <NA>            <NA>         
    ##  8 HVO_RS09905 HVO_1084      TrmB… JK3_cetZ1       <NA>            <NA>         
    ##  9 HVO_RS09910 HVO_1085      TrmB… JK3_cetZ1       <NA>            <NA>         
    ## 10 HVO_RS10145 HVO_1133      TrmB… <NA>            rdfA_ddfA       <NA>         
    ## # ℹ 18 more rows
    ## # ℹ 2 more variables: uniprot <chr>, annotation <chr>

``` r
trmB[trmB$locus_tag %in% tB$locus_tag, ]
```

    ## # A tibble: 26 × 11
    ##    locus_tag   baseMean log2FoldChange lfcSE  stat   pvalue          padj acc   
    ##    <chr>          <dbl>          <dbl> <dbl> <dbl>    <dbl>         <dbl> <chr> 
    ##  1 HVO_RS05875    7290.           2.59 0.254  6.25 4.18e-10 0.0000000485  WP_00…
    ##  2 HVO_RS06270    8517.          -2.20 0.272 -4.41 1.05e- 5 0.000301      WP_00…
    ##  3 HVO_RS06310    3249.           4.07 0.459  6.68 2.35e-11 0.00000000333 WP_00…
    ##  4 HVO_RS06430    6156.          -1.54 0.149 -3.65 2.67e- 4 0.00466       WP_00…
    ##  5 HVO_RS08595   40407.           2.09 0.249  4.40 1.09e- 5 0.000308      WP_00…
    ##  6 HVO_RS09725   32719.           3.30 0.480  4.79 1.65e- 6 0.0000636     WP_00…
    ##  7 HVO_RS09905    8547.           2.65 0.391  4.21 2.59e- 5 0.000644      WP_00…
    ##  8 HVO_RS09910   21357.           2.48 0.342  4.31 1.65e- 5 0.000453      WP_00…
    ##  9 HVO_RS10145   16754.           1.84 0.207  4.04 5.26e- 5 0.00118       WP_00…
    ## 10 HVO_RS10820   26869.           3.31 0.405  5.70 1.20e- 8 0.000000920   WP_00…
    ## # ℹ 16 more rows
    ## # ℹ 3 more variables: old_locus_tag <chr>, annotation <chr>, model <chr>
