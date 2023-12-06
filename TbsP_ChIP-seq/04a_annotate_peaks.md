Annotate peaks
================
Rylee Hackley

Annotating peaks!

``` r
# BiocManager::install(c("GenomicRanges","rtracklayer", "ChIPseeker", "IRanges"))

library(tidyverse)
library(openxlsx)
library(GenomicFeatures)
library(rtracklayer)
```

NC_013964.1 end 437313 NC_013965.1 end 6643 NC_013966.1 start 152
NC_013966.1 end 635564 NC_013967.1 start 258 NC_013967.1 end 2846656
NC_013968.1 start 101 NC_013968.1 end 84525

Load genome file

``` r
gff <- makeTxDbFromGFF("00_genome_files/genomic.gff",
  format = "gff",
  dataSource = "NCBI", organism = "Haloferax volcanii"
)
```

    ## Warning in .extract_transcripts_from_GRanges(tx_IDX, gr, mcols0$type, mcols0$ID, : the transcript names ("tx_name" column in the TxDb object) imported
    ##   from the "Name" attribute are not unique

    ## Warning in makeTxDbFromGRanges(gr, metadata = metadata): The following transcripts were dropped because their exon ranks could
    ##   not be inferred (either because the exons are not on the same
    ##   chromosome/strand or because they are not separated by introns):
    ##   gene-HVO_RS19970

``` r
gff.df <- read_csv("00_genome_files/genomic_gff_key.csv")
```

Get GRanges of genes and promoters and intergenic regions

``` r
genes.only <- genes(gff, columns = c("gene_id"))

# get promoter regions
pro250 <- promoters(genes.only, upstream = 250, downstream = 0) ## Warning: this encroaches on nearby CDS
pro250$type <- rep("promoter", length(pro250))

# merge end of chr info with genes
genes.only$type <- rep("gene", length(genes.only))

# non coding regions of the genome
IG <- gaps(GenomicRanges::reduce(genes.only, ignore.strand = T))
```

overlap promoter region and IG regions to get conservative promoter
ranges (that respect genic regions)

``` r
findOverlaps(pro250, IG, ignore.strand = T) -> pros
pintersect(pro250[queryHits(pros)], IG[subjectHits(pros)]) -> trimmed.pro
trimmed.pro <- trimmed.pro[, 1:2]
names(trimmed.pro) <- NULL

# regenerate IG ranges to respect promoter regions
ig.only <- gaps(GenomicRanges::reduce(c(genes.only, trimmed.pro), ignore.strand = T))

# for motif background models, export bed files for trimmer promoter regions, and all non-coding regions
tmp <- reduce(trimmed.pro, ignore.strand = T)
export.bed(tmp, "05a_motif_discovery/trimmed_promoters.bed")
export.bed(IG, "05a_motif_discovery/intergenic.bed")
```

overlap promoter region and IG regions to get conservative promoter
ranges (that respect genic regions)

``` r
findOverlaps(pro250, IG) -> pros
pintersect(pro250[queryHits(pros)], IG[subjectHits(pros)]) -> trimmed.pro
trimmed.pro <- trimmed.pro[, 1:2]
names(trimmed.pro) <- NULL
```

``` r
(filenames <- list.files(path = "02_DiffBind/", pattern = ".bed"))
```

    ## [1] "02_diffbind_RLE.bed"   "glucose_peaks.bed"     "noglucose_peaks.bed"  
    ## [4] "tbsP_shared_peaks.bed"

``` r
wb <- createWorkbook()

for (i in 1:length(filenames)) {
  # import gff file and convert to dataframe
  bed <- rtracklayer::import(paste("02_DiffBind/", filenames[i], sep = ""), format = "bed", )
  bed.df <- as.data.frame(bed)

  # add informative meta data & order by peak size
  bed.sort <- sort(bed, ~score, decreasing = T, ignore.strand = T) # sort by score
  bed.sort$peakrank <- seq(1, length(bed.sort)) # create peak rank ID by score
  bed.sort$peakwidth <- bed.df$width
  bed.sort[, -1] %>% sort() -> bed # drop names, reorder, but keep peak rank.

  bed.df <- as.data.frame(bed)

  # create hits objects of the overlaps (all peak widths have been standardizes to 300bp wide) adjusting the overlap requirement changes the stringency of my peak annotation.
  GenomicRanges::findOverlaps(genes.only, bed, ignore.strand = T, minoverlap = 1) -> genes
  GenomicRanges::findOverlaps(trimmed.pro, bed, ignore.strand = T, minoverlap = 1) -> promoters
  GenomicRanges::findOverlaps(ig.only, bed, ignore.strand = T, minoverlap = 1) -> ig.regions

  # get IRanges from hits objects and add informative metadata
  genelist <- bed[subjectHits(genes)][, -3]
  genelist$type <- rep("gene", length(genes))
  pintersect(genes.only[queryHits(genes)], bed[subjectHits(genes)]) -> overlaps
  genelist$overlap <- width(overlaps)
  strand(genelist) <- strand(genes.only[queryHits(genes)])
  genelist$TXNAME <- as.character(genes.only$gene_id[queryHits(genes)])

  prolist <- bed[subjectHits(promoters)][, -3]
  prolist$type <- rep("promoter", length(promoters))
  pintersect(trimmed.pro[queryHits(promoters)], bed[subjectHits(promoters)]) -> overlaps
  prolist$overlap <- width(overlaps)
  strand(prolist) <- strand(trimmed.pro[queryHits(promoters)])
  prolist$TXNAME <- as.character(trimmed.pro$gene_id[queryHits(promoters)])

  iglist <- bed[subjectHits(ig.regions)][, -3]
  iglist$type <- rep("intergenic", length(ig.regions))
  pintersect(ig.only[queryHits(ig.regions)], bed[subjectHits(ig.regions)]) -> overlaps
  iglist$overlap <- width(overlaps)
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
  final <- rbind(one, two, three) %>% arrange(seqnames, start, peakrank)
  colnames(final)[c(2, 3, 10)] <- c("peak_start", "peak_end", "locus_tag")

  # merge with gff information (get NCBI annotations and locus names)
  gff.df[gff.df$locus_tag %in% final$locus_tag, ][c(2:7, 10)] -> tmp
  left_join(final, tmp, by = "locus_tag") -> final

  # reorder
  (final <- final[c(7, 1, 2, 3, 4, 6, 8, 9, 10, 12, 5, 14, 15, 13, 16)])
  names(final)[c(1, 2, 6, 7, 8, 11, 12, 13, 14)] <- c("peak_rank", "chr", "peak_score", "overlap_feature", "overlap_length", "feature_strand", "feature_start", "feature_end", "feature_length")
  final$peak_score <- round(final$peak_score, digits = 2)
  final %>% arrange(peak_rank, peak_start) -> final

  # write to excel file, storing results for each .bed as a different sheet
  addWorksheet(wb, str_remove_all(filenames[i], ".bed"), gridLines = TRUE) # Add a worksheet
  writeData(wb, sheet = i, final, rowNames = FALSE) # write data to worksheet i
}
```

    ## Warning in .local(x, decreasing, ...): 'ignore.strand' ignored when 'by' is
    ## specified

    ## Warning in .local(x, decreasing, ...): 'ignore.strand' ignored when 'by' is
    ## specified

    ## Warning in .local(x, decreasing, ...): 'ignore.strand' ignored when 'by' is
    ## specified

    ## Warning in .local(x, decreasing, ...): 'ignore.strand' ignored when 'by' is
    ## specified

``` r
saveWorkbook(wb, "04a_peak_annotation/04a_genelists.xlsx", overwrite = TRUE)
```

in terminal with bedtools loaded, run: bedtools getfasta -fi
../GCF_000223905.1_ASM22390v1_genomic.fna -bed promoters.bed \>
promoter_seqs.fna output beds should be sorted, but if get error, run
bedtools sort -i <BED>

create peak profile

``` r
bed <- rtracklayer::import("03_plot_peaks/ChIP_consensus.bed", format = "bed", )
bed$score <- round(as.numeric(bed$score), digits = 2)
bed$peak_rank <- rank(-bed$score) # peak ID by peak size, congruent with other analysis
bed <- bed[, -1]

GenomicRanges::findOverlaps(genes.only, bed, ignore.strand = T, minoverlap = 0) -> genes
pintersect(bed[subjectHits(genes)], genes.only[queryHits(genes)]) -> g.peaks
g.peaks$hit <- rep("genic", length(g.peaks))
g.peaks$name <- genes.only[queryHits(genes)]$TXNAME
g.peaks$overlap <- width(g.peaks)

GenomicRanges::findOverlaps(trimmed.pro, bed, ignore.strand = T, minoverlap = 0) -> promoters
pintersect(bed[subjectHits(promoters)], trimmed.pro[queryHits(promoters)]) -> p.peaks
p.peaks$hit <- rep("promoter", length(p.peaks))
p.peaks$name <- rep("promoter", length(p.peaks))
p.peaks$overlap <- width(p.peaks)

GenomicRanges::findOverlaps(ig.only, bed, ignore.strand = T, minoverlap = 0) -> intergeneic
pintersect(bed[subjectHits(intergeneic)], ig.only[queryHits(intergeneic)]) -> ig.peaks
ig.peaks$hit <- rep("intergenic", length(ig.peaks))
ig.peaks$name <- rep("intergenic", length(ig.peaks))
ig.peaks$overlap <- width(ig.peaks)

all <- sort(c(g.peaks, p.peaks, ig.peaks))
all <- split(all, all$peak_rank)

# check that each peak exactly the expected width
peak.width <- vector()
for (i in 1:length(all)) {
  tmp <- sum(all[[1]]$overlap)
  peak.width <- append(peak.width, tmp)
}
peak.width # they are! nice!
```

    ## [1] 301 301 301 301 301 301 301 301 301

data wrangle gRangesList to wide df

``` r
as.data.frame(all)[, -c(1:2)] %>%
  pivot_wider(., id_cols = c("peak_rank"), names_from = "hit", values_from = "overlap") -> peaks.final
peaks.final$genic <- sapply(peaks.final$genic, sum)
peaks.final$promoter <- sapply(peaks.final$promoter, sum)
peaks.final$intergenic <- sapply(peaks.final$intergenic, sum)

peaks.final
```

    ## # A tibble: 9 × 4
    ##   peak_rank genic promoter intergenic
    ##       <dbl> <int>    <int>      <int>
    ## 1         1   177      124          0
    ## 2         2   299        2          0
    ## 3         3   111      190          0
    ## 4         4   113      376          0
    ## 5         5   288        0         13
    ## 6         6   125      352          0
    ## 7         7   201      100          0
    ## 8         8     0       26        275
    ## 9         9   102      398          0

``` r
write_csv(peaks.final, "04c_peak_composition/04c_peakcomposition.csv")
write_csv(as.data.frame(all)[, -c(1:2)], "04c_peak_composition/04c_ordered_peakcomp.csv")
```
