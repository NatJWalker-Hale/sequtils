#! /bin/bash
cd "$(dirname "$0")"

# indexes in fasta and maps paired end reads with bwa. Provide, in order, infasta, R1, R2, threads

prefix=$(echo $2 | cut -f1 -d"_")

bwa index $1

bwa mem -t $4 $1 $2 $3 | samtools view -b -f 1 -F 12 > tmp.bam 
samtools fastq -1 "$prefix"_mapped_R1.fq -2 "$prefix"_mapped_R2.fq tmp.bam

gzip "$prefix"_mapped_R1.fq
gzip "$prefix"_mapped_R2.fq

rm tmp.bam
