#!/bin/bash

# extend reads towards their 3´end
for file in *_sorted.bam
do 
  filename=./01_extended_sorted_bam/${file%.*}_extended.bam
  bedtools bamtobed -i $file | bedtools slop -i - -g chrom.sizes.genome -s -r 150 -l 0 | bedtools bedtobam -i - -g chrom.sizes.genome > $filename 
  echo $file extension compled
done

# sort & index
for file in ./01_extended_sorted_bam/*_sorted_extended.bam
do 
  filename=${file/_sorted_extended.bam/}
  filename_out=$filename.sorted.extended.bam
  samtools sort $file -o $filename_out
  samtools index $filename_out
  #echo $file; echo $filename_out
done

#convert bedgraph (for log2 calcs in R)
for file in ./01_extended_sorted_bam/*.sorted.extended.bam
do
  filename=${file/.sorted.extended.bam/}
  name=./02_bedgraph_extended/${filename##*/}.bedgraph
  bedtools genomecov -d -ibam $file > $name
  echo $file position specific genome coverage finsihed
done
