#!/bin/bash

#convert bedgraph (for log2 calcs in R)
for file in *_sorted.bam
do
	name=./00_bedgraph/${file%.*}.bedgraph
	#bedtools genomecov -d -ibam $file > $name
	#echo $file position specific genomecoverage finsihed
        echo $name
done
