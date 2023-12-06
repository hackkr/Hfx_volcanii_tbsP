TbsP shared peaks (trackviewer)
================
Rylee Hackley

``` r
# BiocManager::install(c("Gviz","trackViewer", "AnnotationDbi"))

library(tidyverse)
library(Gviz)
library(GenomicRanges)
library(rtracklayer)
library(GenomicFeatures)
library(AnnotationDbi)
library(trackViewer)
```

Import bam files:

``` r
Filelist <- read_csv("01a_mosaics/01a_hvo_sample_key.csv", col_names = F)

for (i in 1:nrow(Filelist)) {
  assign(paste(Filelist$X3[i], "_bam", sep = ""), importBam(file.path("00_sorted_bams/", Filelist$X1[i])))
}
```

Import bed files of individual reps from mosaics folder:

``` r
filenames <- list.files(path = "01a_mosaics/peaks/", pattern = ".bed")
labels <- str_replace(filenames, ".bed", "_bed")

for (i in 1:length(filenames)) {
  label <- labels[i]
  scores <- importScore(file.path("01a_mosaics/peaks/", filenames[i]), format = "BED")
  assign(label, scores)
}
```

load files for gene model

``` r
gff <- GenomicFeatures::makeTxDbFromGFF("00_genome_files/genomic.gff", format = "gff")
```

    ## Warning in .extract_transcripts_from_GRanges(tx_IDX, gr, mcols0$type, mcols0$ID, : the transcript names ("tx_name" column in the TxDb object) imported
    ##   from the "Name" attribute are not unique

    ## Warning in makeTxDbFromGRanges(gr, metadata = metadata): The following transcripts were dropped because their exon ranks could
    ##   not be inferred (either because the exons are not on the same
    ##   chromosome/strand or because they are not separated by introns):
    ##   gene-HVO_RS19970

``` r
gff.df <- read_csv("00_genome_files/genomic_gff_key.csv")

# create chromosome gRanges
gr <- GRanges("NC_013964.1", IRanges(1, 437313))
gr2 <- GRanges("NC_013966.1", IRanges(1, 635564))
gr3 <- GRanges("NC_013967.1", IRanges(1, 2846656))
gr4 <- GRanges("NC_013968.1", IRanges(1, 84525))
```

Sizes of HVO genomic elements: NC_013964.1 end 437313 NC_013965.1 end
6643 NC_013966.1 end 635564 NC_013967.1 end 2846656 NC_013968.1 end
84525

build gene model

``` r
# create bed files of genes with strand info
some.genes <- GenomicFeatures::transcripts(gff)
hvo.genes.plus <- some.genes[strand(some.genes) == "+"]
hvo.genes.plus$score <- rep(1, length(hvo.genes.plus))

hvo.genes.minus <- some.genes[strand(some.genes) == "-"]
hvo.genes.minus$score <- rep(1, length(hvo.genes.minus))

# export beds to import later...
rtracklayer::export.bed(hvo.genes.minus, "03_plot_peaks/hvo_genes_minus.bed")
rtracklayer::export.bed(hvo.genes.plus, "03_plot_peaks/hvo_genes_plus.bed")

# import both files to make stranded gene and motif tracks
hvo.genes.bed <- importScore(file.path("03_plot_peaks/hvo_genes_plus.bed"), file.path("03_plot_peaks/hvo_genes_minus.bed"), format = "BED")
strand(hvo.genes.bed$dat) <- "+"
strand(hvo.genes.bed$dat2) <- "-"

# asRNAs from Gelsinger 2018
as.rna <- importScore(file.path("03_plot_peaks/asRNA_plus.bed"),
  file.path("03_plot_peaks/asRNA_minus.bed"),
  format = "BED"
)
strand(as.rna$dat) <- "+"
strand(as.rna$dat2) <- "-"
```

import beds of peaks shared across bioreps from 02_peaklists folder:

``` r
average.scores <- read_delim("02_DiffBind/tbsP_shared_peaks.bed",
  delim = "\t",
  col_names = c("chr", "start", "end", "strand", "score")
)[, 5]
max(average.scores$score)
```

    ## [1] 821.6865

``` r
hvo.tbsP <- importScore("02_DiffBind/tbsP_shared_peaks.bed", format = "BED")
hvo.tbsP$dat <- coverageGR(hvo.tbsP$dat)
```

``` r
gr.goi <- rtracklayer::import("02_DiffBind/tbsP_shared_peaks.bed", format = "bed")

# graph each GOI individually and save in directory by category
for (i in 1:length(gr.goi)) {
  filename <- paste("03_plot_peaks/trackviewer_shared/", "peak_num_",
    gr.goi$name[i], ".png",
    sep = ""
  )
  title <- paste("peak_num_", gr.goi$name[i], sep = "")
  begin <- start(gr.goi[i])
  term <- end(gr.goi[i])

  # set figure range around peak
  if ((begin - 1000) < 1) {
    term2 <- term + 1000
    gr.tmp <- GRanges(as.vector(seqnames(gr.goi[i])), IRanges(begin, term2))
  } else {
    gr.tmp <- GRanges(as.vector(seqnames(gr.goi[i])), IRanges(begin, term)) + 1000
  }

  # create list of desired tracks and name
  trackList <- trackList(as.rna, hvo.genes.bed, hvo.tbsP, HVO_tbsPHA_2_bam, HVO_tbsPHA_glc_2_bam)
  names(trackList) <- c("anti-sense RNA", "genes", "consensus peak", "no glucose", "glucose")

  optSty <- optimizeStyle(trackList, theme = "safe")
  trackList <- optSty$tracks
  viewerStyle <- optSty$style

  # set margins and individual track heights. Tracks are plotted from bottom up.
  setTrackViewerStyleParam(viewerStyle, "margin", c(.06, .09, .01, .09))
  setTrackStyleParam(trackList[[5]], "height", .35)
  setTrackStyleParam(trackList[[4]], "height", .35)
  setTrackStyleParam(trackList[[3]], "height", .1)
  setTrackStyleParam(trackList[[2]], "height", .1)
  setTrackStyleParam(trackList[[1]], "height", .1)

  setTrackStyleParam(trackList[[3]], "ylim", c(0, 1))
  setTrackStyleParam(trackList[[2]], "ylim", c(-1, 1))
  setTrackStyleParam(trackList[[1]], "ylim", c(-1, 1))

  # set track label positions and sizes
  for (x in 1:length(trackList)) {
    setTrackStyleParam(trackList[[x]], "ylabgp", list(cex = 1))
    setTrackStyleParam(trackList[[x]], "ylabpos", "topleft")
    setTrackStyleParam(trackList[[x]], "marginTop", .2)
  }

  # set track colors
  setTrackStyleParam(trackList[[5]], "color", "skyblue3")
  setTrackStyleParam(trackList[[4]], "color", "grey25")
  setTrackStyleParam(trackList[[3]], "color", "grey")
  setTrackStyleParam(trackList[[2]], "color", c("light blue", "light blue"))
  setTrackStyleParam(trackList[[1]], "color", c("mediumpurple", "mediumpurple"))

  # save plots with title and subtitle
  png(filename = filename, width = 800, height = 1000)
  viewTracks(trackList, gr = gr.tmp, viewerStyle = viewerStyle)
  dev.off()
}
```
