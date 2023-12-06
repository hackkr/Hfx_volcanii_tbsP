05_annotate_FIMO
================
Rylee Hackley

``` r
library(tidyverse)
library(readxl)
library(GenomicFeatures)
library(janitor)
```

``` r
motifs <- read_csv("05_promoter_MEME/trimmed_1st/1st_FIMO.csv")
gr <- makeGRangesFromDataFrame(motifs, seqnames.field = "sequence_name", keep.extra.columns = T)

gff <- makeTxDbFromGFF("000_20230911_hvo/genomic.gff", format = "gff", dataSource = "NCBI", organism = "Haloferax volcanii")
gff.df <- read_csv("000_20230911_hvo/genomic_gff_key.csv")
hvo.cogs <- read.delim("000_20230911_hvo/arcogs-14-18.hvo.updated.schronheit.txt", sep = "\t")
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
IG <- gaps(IRanges::reduce(genes.only, ignore.strand = T))
```

overlap promoter region and IG regions to get conservative promoter
ranges (that respect genic regions)

``` r
findOverlaps(pro250, IG) -> pros
pintersect(pro250[queryHits(pros)], IG[subjectHits(pros)]) -> trimmed.pro
trimmed.pro <- trimmed.pro[, 1:2]
names(trimmed.pro) <- NULL
```

regenerate IG ranges to respect promoter regions

``` r
ig.only <- gaps(IRanges::reduce(c(genes.only, trimmed.pro), ignore.strand = T))
```

\#Annotate FIMO results: “where are the motif hits across the
genome?” 1. identify and annotate all motifs that are in
intergenic/promoter regions (use p-value/score cutoff)  
max(peak.with.motif\$motif.score, na.rm = T) 3. look for pairs/groups of
motifs that overlap -\> there are 3, and included in the results file
for all fimo hits

``` r
GenomicRanges::findOverlaps(genes.only, gr, ignore.strand = T, minoverlap = 5) -> genes
GenomicRanges::findOverlaps(trimmed.pro, gr, ignore.strand = T, minoverlap = 5) -> promoters
GenomicRanges::findOverlaps(ig.only, gr, ignore.strand = T, minoverlap = 5) -> ig.regions

# get IRanges from hits objects and add informative metadata
genelist <- gr[subjectHits(genes)][, -3]
genelist$type <- rep("gene", length(genes))
strand(genelist) <- strand(genes.only[queryHits(genes)])
genelist$feature_start <- start(genes.only[queryHits(genes)])
genelist$feature_end <- end(genes.only[queryHits(genes)])
genelist$TXNAME <- as.character(genes.only$GENEID[queryHits(genes)])

prolist <- gr[subjectHits(promoters)][, -3]
prolist$type <- rep("promoter", length(promoters))
prolist$feature_start <- start(trimmed.pro[queryHits(promoters)])
prolist$feature_end <- end(trimmed.pro[queryHits(promoters)])
strand(prolist) <- strand(trimmed.pro[queryHits(promoters)])
prolist$TXNAME <- as.character(trimmed.pro$GENEID[queryHits(promoters)])

iglist <- gr[subjectHits(ig.regions)][, -3]
iglist$type <- rep("intergenic", length(ig.regions))
iglist$feature_start <- start(ig.only[queryHits(ig.regions)])
iglist$feature_end <- end(ig.only[queryHits(ig.regions)])
iglist$TXNAME <- NA

# convert separate IRanges to Data frames
seqs <- seq(1, length(genes))
as.data.frame(prolist) -> one
rownames(one) <- NULL
as.data.frame(genelist, row.names(seqs)) -> two
rownames(two) <- NULL
as.data.frame(iglist, row.names(seqs)) -> three
rownames(three) <- NULL

# combine dfs (gene hits and promoter hits)
final <- rbind(one, two, three) %>% distinct(.keep_all = T)
colnames(final)[c(2, 3, 14)] <- c("motif_start", "motif_end", "locus_tag")

# merge with gff information (get NCBI annotations and locus names)
gff.df[gff.df$locus_tag %in% final$locus_tag, ] -> tmp
tmp[c(2, 3, 4, 10)] -> tmp2
left_join(final, tmp2, by = "locus_tag") %>% arrange(p.value) -> final

# check that all are accounted for:
nrow(final) == nrow(one) + nrow(two) + nrow(three)
```

    ## [1] TRUE

``` r
# now do p-value cutoff
final %>% filter(q.value < 1) -> final

write_csv(final, "05_promoter_MEME/1st_annotated_FIMO_all.csv")
```

``` r
final$p.value %>% max()
```

    ## [1] 5.73e-06

``` r
final$q.value %>% max()
```

    ## [1] 0.958

``` r
final$locus_tag %>%
  unique() %>%
  length()
```

    ## [1] 52

``` r
filter(final, final$type == "promoter") %>%
  unique() %>%
  length()
```

    ## [1] 17
