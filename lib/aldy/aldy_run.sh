#!/bin/bash

gene=$1
sampleNo=$2
bamPath=$3
outPath=$4

#for i in $(seq -f "%03g" 1 48)
#do
mkdir -p $outPath/$gene
aldy genotype -p illumina -g $gene -o $outPath/$gene/$sampleNo.$gene.aldy $bamPath --log $outPath/$gene/aldy.log #$sampleNo.realigned.bam
#done
