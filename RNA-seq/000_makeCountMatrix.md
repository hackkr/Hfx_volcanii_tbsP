makeCountMatrix
================
Rylee Hackley
OCT 2023

``` r
library(tidyverse)
library(magrittr)
library(purrr)
library(DESeq2)
```

``` r
#read all counts (from 2019) into single dataframe
counts2019 <- read_tsv("00_counts/2019_counts.txt", skip = 1)
tmp <- colnames(counts2019)
tmp2 <- str_replace_all(tmp, "-", "_")
tmp3 <- str_remove(tmp2, "_sorted.bam")
colnames(counts2019) <- tmp3

## remove counts with less than 60% alignment
tmp <- c("dtbsP_4hmfe_5_S46", "dtbsP_4hmfe_2_S40", "dtrmBdtbsP_4hmfe_2_S56", "dtrmBdtbsP_4hpfe_2_S54", "dtrmBdtbsP_4hpfe_6_S69")

counts2019 %<>% select(!tmp)

#check genotype matches label:
counts2019 %>%
  filter(Geneid == "HVO_RS18535" | Geneid == "HVO_RS17680")
```

    ## # A tibble: 2 × 25
    ##   Geneid   Chr    Start    End Strand Length dpyrE2_4hmfe_2_S5 dpyrE2_4hmfe_4_S9
    ##   <chr>    <chr>  <dbl>  <dbl> <chr>   <dbl>             <dbl>             <dbl>
    ## 1 HVO_RS1… NC_0… 2.53e6 2.54e6 +        1062               301               907
    ## 2 HVO_RS1… NC_0… 2.70e6 2.70e6 -         828                16              2129
    ## # ℹ 17 more variables: dpyrE2_4hmfe_5_S18 <dbl>, dpyrE2_4hpfe_2_S2 <dbl>,
    ## #   dpyrE2_4hpfe_4_S16 <dbl>, dpyrE2_4hpfe_5_S12 <dbl>,
    ## #   dtbsP_4hmfe_6_S52 <dbl>, dtbsP_4hpfe_2_S38 <dbl>, dtbsP_4hpfe_5_S43 <dbl>,
    ## #   dtbsP_4hpfe_6_S49 <dbl>, dtrmB_4hmfe_1_S24 <dbl>, dtrmB_4hmfe_5_S30 <dbl>,
    ## #   dtrmB_4hmfe_6_S36 <dbl>, dtrmB_4hpfe_1_S21 <dbl>, dtrmB_4hpfe_5_S27 <dbl>,
    ## #   dtrmB_4hpfe_6_S33 <dbl>, dtrmBdtbsP_4hmfe_5_S66 <dbl>,
    ## #   dtrmBdtbsP_4hmfe_6_S72 <dbl>, dtrmBdtbsP_4hpfe_5_S65 <dbl>

``` r
#remove mislabeled samples
tmp <- c("dtbsP_4hmfe_6_S52", "dtrmBdtbsP_4hpfe_5_S65") #(have tbsP and trmB expr when shouldn't)
counts2019 %<>% select(!tmp)

counts2019wt <- read_tsv("00_counts/2019_counts2.txt", skip = 1)[16]
colnames(counts2019wt) <- str_replace_all(colnames(counts2019wt), "-", "_")
colnames(counts2019wt) <- str_remove(colnames(counts2019wt), "_L002_R1_001.fastq.gz_sorted.bam")

counts2019 <- cbind(counts2019, counts2019wt)
```

``` r
counts2021 <- read_tsv("00_counts/2021_counts.txt", skip = 1)

#drop mislabeled hvo sample in column 15:
counts2021 <- counts2021[-15]
colnames(counts2021)[7] <- "HVO_DS2_reference"

colnames(counts2021) <- str_remove(colnames(counts2021) , "_sorted.bam")

#check genotype matches label:
counts2021 %>%
  filter(Geneid == "HVO_RS18535" | Geneid == "HVO_RS17680")
```

    ## # A tibble: 2 × 37
    ##   Geneid      Chr           Start     End Strand Length HVO_DS2_reference
    ##   <chr>       <chr>         <dbl>   <dbl> <chr>   <dbl>             <dbl>
    ## 1 HVO_RS17680 NC_013967.1 2534459 2535520 +        1062              4316
    ## 2 HVO_RS18535 NC_013967.1 2699406 2700233 -         828              9745
    ## # ℹ 30 more variables: deltadelta_a_glu_S16 <dbl>,
    ## #   deltadelta_a_noglu_S35 <dbl>, deltadelta_b_glu_S53 <dbl>,
    ## #   deltadelta_b_noglu_S36 <dbl>, deltadelta_c_glu_S43 <dbl>,
    ## #   deltadelta_d_glu_S44 <dbl>, deltadelta_d_noglu_S5 <dbl>,
    ## #   tbsP_a_glu_S33 <dbl>, tbsP_a_noglu_S24 <dbl>, tbsP_b_glu_S2 <dbl>,
    ## #   tbsP_b_noglu_S30 <dbl>, tbsP_c_glu_S51 <dbl>, tbsP_c_noglu_S15 <dbl>,
    ## #   tbsP_d_glu_S52 <dbl>, tbsP_d_noglu_S31 <dbl>, tbsP_e_glu_S3 <dbl>, …

``` r
#correct swapped trmB_a_glu and tbsP_e_glu samples: 
a <- colnames(counts2021)[23]
b <- colnames(counts2021)[25]

colnames(counts2021)[23] <- b
colnames(counts2021)[25] <- a

#merge
(left_join(counts2021, counts2019, by = c("Geneid", "Chr", "Start", "End", "Strand", "Length")) -> counts_combined)
```

    ## # A tibble: 3,982 × 55
    ##    Geneid Chr   Start   End Strand Length HVO_DS2_reference deltadelta_a_glu_S16
    ##    <chr>  <chr> <dbl> <dbl> <chr>   <dbl>             <dbl>                <dbl>
    ##  1 HVO_R… NC_0…   258  1952 +        1695              6411                12475
    ##  2 HVO_R… NC_0…  2000  2797 -         798              8166                12123
    ##  3 HVO_R… NC_0…  2942  4543 +        1602              7809                22141
    ##  4 HVO_R… NC_0…  4686  5813 +        1128              2023                 3800
    ##  5 HVO_R… NC_0…  5936  6745 -         810              2168                 2848
    ##  6 HVO_R… NC_0…  6854  7282 -         429              1087                 1715
    ##  7 HVO_R… NC_0…  7470  8648 +        1179              4786                 8816
    ##  8 HVO_R… NC_0…  8873 10219 +        1347              7631                 1733
    ##  9 HVO_R… NC_0… 10248 10907 -         660              2042                 2390
    ## 10 HVO_R… NC_0… 10956 11528 -         573               747                 4303
    ## # ℹ 3,972 more rows
    ## # ℹ 47 more variables: deltadelta_a_noglu_S35 <dbl>,
    ## #   deltadelta_b_glu_S53 <dbl>, deltadelta_b_noglu_S36 <dbl>,
    ## #   deltadelta_c_glu_S43 <dbl>, deltadelta_d_glu_S44 <dbl>,
    ## #   deltadelta_d_noglu_S5 <dbl>, tbsP_a_glu_S33 <dbl>, tbsP_a_noglu_S24 <dbl>,
    ## #   tbsP_b_glu_S2 <dbl>, tbsP_b_noglu_S30 <dbl>, tbsP_c_glu_S51 <dbl>,
    ## #   tbsP_c_noglu_S15 <dbl>, tbsP_d_glu_S52 <dbl>, tbsP_d_noglu_S31 <dbl>, …

``` r
#calc total reads per sample
as.data.frame(colSums(counts_combined[-c(1:6)])) -> cnt_sums
cnt_sums <- rownames_to_column(cnt_sums, "sample_name")

#export combined counts
write_csv(counts_combined, "00_counts/00_combined_counts.csv")
write_csv(counts_combined, "00_deseq2_input/00_combined_data.csv")
write_csv(cnt_sums, "00_counts/00_counts_sums.csv")
```
