#!/bin/sh

rm data/*

seqretsplit -sequence $1 -outseq .fasta
mv *.fasta data/

for file in data/*.fasta; do perl src/DiTriNucFreq.pl $file $2; done

rm data/*.temp
