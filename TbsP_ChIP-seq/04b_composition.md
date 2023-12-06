binomial test peak location
================
Rylee Hackley

Annotating peaks

``` r
# BiocManager::install(c("GenomicRanges","rtracklayer", "ChIPseeker", "IRanges"))

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
library(openxlsx)
library(GenomicFeatures)
```

    ## Loading required package: BiocGenerics
    ## 
    ## Attaching package: 'BiocGenerics'
    ## 
    ## The following objects are masked from 'package:lubridate':
    ## 
    ##     intersect, setdiff, union
    ## 
    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     combine, intersect, setdiff, union
    ## 
    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs
    ## 
    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    ##     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
    ##     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
    ##     Position, rank, rbind, Reduce, rownames, sapply, setdiff, sort,
    ##     table, tapply, union, unique, unsplit, which.max, which.min
    ## 
    ## Loading required package: S4Vectors
    ## Loading required package: stats4
    ## 
    ## Attaching package: 'S4Vectors'
    ## 
    ## The following objects are masked from 'package:lubridate':
    ## 
    ##     second, second<-
    ## 
    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     first, rename
    ## 
    ## The following object is masked from 'package:tidyr':
    ## 
    ##     expand
    ## 
    ## The following object is masked from 'package:utils':
    ## 
    ##     findMatches
    ## 
    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname
    ## 
    ## Loading required package: IRanges
    ## 
    ## Attaching package: 'IRanges'
    ## 
    ## The following object is masked from 'package:lubridate':
    ## 
    ##     %within%
    ## 
    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     collapse, desc, slice
    ## 
    ## The following object is masked from 'package:purrr':
    ## 
    ##     reduce
    ## 
    ## The following object is masked from 'package:grDevices':
    ## 
    ##     windows
    ## 
    ## Loading required package: GenomeInfoDb
    ## Loading required package: GenomicRanges
    ## Loading required package: AnnotationDbi
    ## Loading required package: Biobase
    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.
    ## 
    ## 
    ## Attaching package: 'AnnotationDbi'
    ## 
    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select

``` r
library(rtracklayer)
```

Load genome file

``` r
gff <- makeTxDbFromGFF("00_genome_files/genomic.gff",
  format = "gff",
  dataSource = "NCBI", organism = "Haloferax volcanii"
)
gff_df <- read_csv("00_genome_files/genomic_gff_key.csv")
consensus <- read_csv("03_plot_peaks/ChIP_consensus.csv")
consensus_gr <- makeGRangesFromDataFrame(consensus, keep.extra.columns = T)
```

Get GRanges of genes and promoters and intergenic regions

``` r
genes.only <- genes(gff, columns = c("GENEID"))

# get promoter regions
pro250 <- promoters(genes.only, upstream = 250, downstream = 0) ## Warning: this encroaches on nearby CDS
pro250$type <- rep("promoter", length(pro250))

# merge end of chr info with genes
genes.only$type <- rep("gene", length(genes.only))

# non coding regions of the genome
IG <- gaps(GenomicRanges::reduce(genes.only, ignore.strand = T))

# overlap promoter region and IG regions to get conservative promoter ranges (that respect genic regions)
findOverlaps(pro250, IG, ignore.strand = T) -> pros
pintersect(pro250[queryHits(pros)], IG[subjectHits(pros)]) -> trimmed.pro
trimmed.pro <- trimmed.pro[, 1:2]
names(trimmed.pro) <- NULL

# regenerate IG ranges to respect promoter regions
ig.only <- gaps(GenomicRanges::reduce(c(genes.only, trimmed.pro), ignore.strand = T))
```

create peak profile

``` r
GenomicRanges::findOverlaps(genes.only, consensus_gr, minoverlap = 1) -> genes
pintersect(consensus_gr[subjectHits(genes)], genes.only[queryHits(genes)]) -> g.peaks
g.peaks$hit <- rep("genic", length(g.peaks))
g.peaks$name <- as.character(genes.only[queryHits(genes)]$GENEID)

GenomicRanges::findOverlaps(trimmed.pro, consensus_gr, ignore.strand = T, minoverlap = 1) -> promoters
pintersect(consensus_gr[subjectHits(promoters)], trimmed.pro[queryHits(promoters)]) -> p.peaks
GenomicRanges::reduce(p.peaks, with.revmap = T, ignore.strand = T) -> p.peaks2
p.peaks2$score <- p.peaks$score[sapply(p.peaks2$revmap, "[", 1)]
p.peaks2$hit <- rep("promoter", length(p.peaks2))
p.peaks2$name <- rep("promoter", length(p.peaks2))
p.peaks <- p.peaks2[, -1]

GenomicRanges::findOverlaps(ig.only, consensus_gr, ignore.strand = T, minoverlap = 1) -> intergeneic
pintersect(consensus_gr[subjectHits(intergeneic)], ig.only[queryHits(intergeneic)]) -> ig.peaks
ig.peaks$hit <- rep("intergenic", length(ig.peaks))
ig.peaks$name <- rep("intergenic", length(ig.peaks))

all <- sort(c(g.peaks, p.peaks, ig.peaks), ignore.strand = T)
all <- split(all, all$score)

# check that each peak exactly the expected width
peak.width <- vector()
for (i in 1:length(all)) {
  tmp <- sum(width(all[[i]]))
  peak.width <- append(peak.width, tmp)
}
peak.width
```

    ## [1] 302 302 302 302 302 302 302 302 302

data wrangle gRangesList to wide df

``` r
as.data.frame(all)[, -c(1:2)] %>%
  pivot_wider(., id_cols = "score", names_from = "hit", values_from = "width") -> peaks.final
```

    ## Warning: Values from `width` are not uniquely identified; output will contain list-cols.
    ## • Use `values_fn = list` to suppress this warning.
    ## • Use `values_fn = {summary_fun}` to summarise duplicates.
    ## • Use the following dplyr code to identify duplicates.
    ##   {data} %>%
    ##   dplyr::group_by(score, hit) %>%
    ##   dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    ##   dplyr::filter(n > 1L)

``` r
peaks.final$genic <- sapply(peaks.final$genic, sum)
peaks.final$promoter <- sapply(peaks.final$promoter, sum)
peaks.final$intergenic <- sapply(peaks.final$intergenic, sum)
# peaks.final$motif <- sapply(peaks.final$motif, sum)

all.df <- as.data.frame(all)[, -c(1:2)]
names(all.df$name) <- NULL
all.df$name <- unlist(all.df$name)

write_csv(all.df, "04c_peak_composition/04c_ordered_peak_composition_consensus.csv")
```

NC_013967.1 2847757 NC_013968.1 85092 NC_013965.1 6359 NC_013964.1
437906 NC_013966.1 635786

``` r
genome.size <- sum(2847757, 85092, 437906, 635786)
sum.all.genes <- sum(IRanges::width(reduce(genes.only)))

all.df %>%
  filter(hit == "intergenic") %>%
  tally(width) -> sum.ig.pks
all.df %>%
  filter(hit == "genic") %>%
  tally(width) -> sum.gene.peaks
all.df %>%
  filter(hit == "promoter") %>%
  tally(width) -> sum.pro.pks

# get expected percentage of reads in genes based on gene density of genome
out <- data.frame(
  "no.peaks" = length(unique(all.df$score)),
  "bp in peaks" = sum(all.df$width),
  "expected in" = sum.all.genes / genome.size, # expected proportion of peaks in genes based on proportion of genes in genome
  "expected out" = 1 - (sum.all.genes / genome.size),
  "obs. in" = (sum.gene.peaks$n) / sum(all.df$width),
  "obs. out" = (sum.pro.pks$n + sum.ig.pks$n) / sum(all.df$width)
)
out
```

    ##   no.peaks bp.in.peaks expected.in expected.out   obs..in  obs..out
    ## 1        9        2718    0.866251     0.133749 0.5235467 0.4764533

``` r
# expected 86.6% of bases in peaks to overlap genes, observed 52%
# expected 13.3% of bases in peaks to overlap intergenic regions, observed 48%

# binomial test, is the proportion of genic bps in peaks greater than we'd expect?
binom.test(sum.gene.peaks$n, out$bp.in.peaks, out$expected.in, alternative = "greater") # no
```

    ## 
    ##  Exact binomial test
    ## 
    ## data:  sum.gene.peaks$n and out$bp.in.peaks
    ## number of successes = 1423, number of trials = 2718, p-value = 1
    ## alternative hypothesis: true probability of success is greater than 0.866251
    ## 95 percent confidence interval:
    ##  0.5075919 1.0000000
    ## sample estimates:
    ## probability of success 
    ##              0.5235467

``` r
binom.test(sum.gene.peaks$n, out$bp.in.peaks, out$expected.in, alternative = "less") # its much less than we'd expect
```

    ## 
    ##  Exact binomial test
    ## 
    ## data:  sum.gene.peaks$n and out$bp.in.peaks
    ## number of successes = 1423, number of trials = 2718, p-value < 2.2e-16
    ## alternative hypothesis: true probability of success is less than 0.866251
    ## 95 percent confidence interval:
    ##  0.000000 0.539464
    ## sample estimates:
    ## probability of success 
    ##              0.5235467

``` r
# binomial test, is the proportion of bps in peaks falling outside of gene sequences greater than we'd expect?
binom.test(sum.ig.pks$n + sum.pro.pks$n, out$bp.in.peaks, out$expected.out, alternative = "greater") # yes
```

    ## 
    ##  Exact binomial test
    ## 
    ## data:  sum.ig.pks$n + sum.pro.pks$n and out$bp.in.peaks
    ## number of successes = 1295, number of trials = 2718, p-value < 2.2e-16
    ## alternative hypothesis: true probability of success is greater than 0.133749
    ## 95 percent confidence interval:
    ##  0.460536 1.000000
    ## sample estimates:
    ## probability of success 
    ##              0.4764533
