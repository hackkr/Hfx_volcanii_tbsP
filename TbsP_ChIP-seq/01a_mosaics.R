### Peak calling script using MOSAiCS R package
# Cynthia L. Darnell and Amy K. Schmid

# MOSAiCS 
# Dongjun Chung, Pei Fen Kuan, Rene Welch, Sunduz Keles
# https://bioconductor.org/packages/release/bioc/html/mosaics.html

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("mosaics")

library(mosaics)
library(hexbin)
library(tidyverse)

sample_file <- read_csv("01a_mosaics/01a_hvo_sample_key.csv", col_names = F)
IP_files <- unique(sample_file$X1)
WCE_files <- unique(sample_file$X2)


#construct bins
for (i in 1:nrow(sample_file)){
  constructBins(infile=paste("00_sorted_bams/", IP_files[i], sep=""),
                fileFormat="bam", outfileLoc="01a_mosaics/bins/", byChr=FALSE,
                fragLen=200, binSize=200, capping=0, PET=FALSE)
  
  constructBins(infile=paste("00_sorted_bams/", WCE_files[i], sep=""),
                fileFormat="bam", outfileLoc="01a_mosaics/bins/", byChr=FALSE,
                fragLen=200, binSize=200, capping=0, PET=FALSE)
}

for (i in 1:nrow(sample_file)) {
  sample_name <- paste("01a_mosaics/bins/", sample_file[i,1], sep = "")
  sample_name <- str_replace(string = sample_name, pattern = ".bam", replacement = ".bam_fragL200_bin200.txt")
  ref_name <- paste("01a_mosaics/bins/", sample_file[i,2], sep = "")
  ref_name <- str_replace(string = ref_name, pattern = ".bam", replacement = ".bam_fragL200_bin200.txt")
  
  print(paste("analyzing", sample_name, "against", ref_name))
  
  binTest <- readBins(type=c("chip", "input"), fileName= c(sample_name, ref_name))
  plot(binTest)
  plot(binTest, plotType="input")
  dev.copy(png, paste("01a_mosaics/peaks/plots/", sample_file$X3[i], "_tagcounts.png", sep=""))
  dev.off()
  
  count_data <- hexbin (binTest@input, binTest@tagCount, xbins=100)
  control <- plot(count_data, trans=log, inv=exp, colramp=rainbow, xlab="WCE", ylab="ChIP", lcex=0.9)
  hexVP.abline(control$plot.vp, a=0, b=sum(binTest@tagCount)/sum(binTest@input), lwd=0.2)
  dev.copy(png, paste("01a_mosaics/peaks/plots/", sample_file$X3[i], "_counts.png", sep=""))
  dev.off()
  
  fitTest <- mosaicsFit(binTest, analysisType="IO", bgEst="automatic")
  plot(fitTest)
  dev.copy(png, paste("01a_mosaics/peaks/plots/", sample_file$X3[i], "_fit.png", sep=""))
  dev.off()
  
  peakTest <- mosaicsPeak(fitTest, signalModel="2S", FDR=0.05)
  export(peakTest, type="bed", filename=paste("01a_mosaics/peaks/", sample_file$X3[i], ".bed", sep=""))
  export(peakTest, type="txt", filename=paste("01a_mosaics/peaks/", sample_file$X3[i],  ".txt", sep=""))
}

#Initial nearby peaks are merged if the distance (in bp) between them is less than 'maxgap'. 
#Some initial peaks are removed if their lengths are shorter than 'minsize'.

#If you use a bin size shorter than the average fragment length in the experiment, 
  #we recommend to set 'maxgap' to the average fragment length and 'minsize' to the bin size. 
  #This setting removes peaks that are too narrow (e.g., singletons). 
#If you set the bin size to the average fragment length (or maybe bin size is larger than the average fragment length), 
  #we recommend setting 'minsize' to a value smaller than the average fragment length while leaving 'maxgap' the same as the average fragment length. 
  #This is to prevent filtering using 'minsize' because initial peaks would already be at a reasonable width. 
