02_schiller_gene_list
================
Rylee Hackley

``` r
library(tidyverse)
```

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
library(readxl)

supp2 <- read_xlsx("02_schiller_overlap/2023_biorxiv_supp/media-2.xlsx", sheet = 3) # sig abuncance difference for wild type between early log and late log as well as at least one shape-specific comparison
supp3 <- read_xlsx("02_schiller_overlap/2023_biorxiv_supp/media-3.xlsx")[c(2, 9)] %>% # clusters
  filter(cluster == 11 | cluster == 12 | cluster == 14)
```

    ## Warning: Coercing text to numeric in H2097 / R2097C8: 'NaN'

    ## New names:
    ## • `` -> `...1`

``` r
supp4 <- read_xlsx("02_schiller_overlap/2023_biorxiv_supp/media-4.xlsx", sheet = 3) # sig abundance differences for wild type between early log and late log as well as at least one shape-specific comparison

annotations <- read_xlsx("02_schiller_overlap/2023_biorxiv_supp/media-3.xlsx")[9]
```

    ## Warning: Coercing text to numeric in H2097 / R2097C8: 'NaN'

    ## New names:
    ## • `` -> `...1`

wrangle to get locus tags

``` r
supp2 %>%
  separate_longer_delim(`Protein ID`, delim = "<|>") %>%
  separate_wider_delim(`Protein ID`, delim = " [", names = c("old_locus_tag", "annotation")) %>%
  separate_wider_delim(annotation, delim = "] ", names = c("uniprot", "annotation")) %>%
  mutate(analysis_type = "JK3_cetZ1") -> supp2.1

supp3 %>%
  separate_longer_delim("names", delim = "<|>") %>%
  separate_wider_delim(names, delim = " [", names = c("old_locus_tag", "annotation")) %>%
  separate_wider_delim(annotation, delim = "] ", names = c("uniprot", "annotation")) %>%
  mutate(analysis_type = paste("clusters", cluster, sep = "_")) -> supp3.1

annotations %>%
  separate_longer_delim("names", delim = "<|>") %>%
  separate_wider_delim(names, delim = " [", names = c("old_locus_tag", "annotation"), too_many = "merge") %>%
  separate_wider_delim(annotation, delim = "] ", names = c("uniprot", "annotation"), too_many = "merge") %>%
  distinct() -> annotations1

supp4 %>%
  separate_longer_delim(`Protein ID`, delim = "<|>") %>%
  separate_wider_delim(`Protein ID`, delim = " [", names = c("old_locus_tag", "annotation")) %>%
  separate_wider_delim(annotation, delim = "] ", names = c("uniprot", "annotation")) %>%
  mutate(analysis_type = "rdfA_ddfA") -> supp4.1
```

``` r
full_join(supp2.1[c(1:3, 20)], supp4.1[c(1:3, 14)], by = c("old_locus_tag", "uniprot", "annotation")) -> tmp

full_join(tmp, supp3.1[c(2:5)], by = c("old_locus_tag", "uniprot", "annotation")) -> all.genes
```

    ## Warning in full_join(tmp, supp3.1[c(2:5)], by = c("old_locus_tag", "uniprot", : Detected an unexpected many-to-many relationship between `x` and `y`.
    ## ℹ Row 35 of `x` matches multiple rows in `y`.
    ## ℹ Row 73 of `y` matches multiple rows in `x`.
    ## ℹ If a many-to-many relationship is expected, set `relationship =
    ##   "many-to-many"` to silence this warning.

``` r
all.genes[c(1, 4:6)] %>%
  pivot_longer(!old_locus_tag, names_to = "supp_file", values_to = "analysis") %>%
  drop_na() -> all.genes2

all.genes2 %>%
  pivot_wider(names_from = supp_file, values_from = analysis) %>%
  unnest() %>%
  distinct() %>%
  arrange(old_locus_tag) -> all.genes2
```

    ## Warning: Values from `analysis` are not uniquely identified; output will contain
    ## list-cols.
    ## • Use `values_fn = list` to suppress this warning.
    ## • Use `values_fn = {summary_fun}` to summarise duplicates.
    ## • Use the following dplyr code to identify duplicates.
    ##   {data} %>%
    ##   dplyr::group_by(old_locus_tag, supp_file) %>%
    ##   dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    ##   dplyr::filter(n > 1L)

    ## Warning: `cols` is now required when using `unnest()`.
    ## ℹ Please use `cols = c(analysis_type.x, analysis_type.y, analysis_type)`.

``` r
all.genes3 <- inner_join(all.genes2, annotations1)
```

    ## Joining with `by = join_by(old_locus_tag)`

``` r
write_csv(all.genes3, "02_schiller_overlap/02_schiller2023_list.csv")
```
