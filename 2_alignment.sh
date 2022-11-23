#!/bin/bash

#SBATCH -A e30836
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -t 12:00:00
#SBATCH --mem=80G

module load STAR/2.5.2

input="./srr-nums"

while read -r line
do
	mkdir /projects/e30836/project/group2-f22/alignment/$line
	cd /projects/e30836/project/group2-f22/alignment/$line
	
	STAR --runThreadN 10 --quantMode TranscriptomeSAM --genomeDir /projects/e30836/hw1/hg38.index/STAR --readFilesIn /projects/e30836/project/group2-f22/fastq-data/$line/"$line"_1.fastq
done < "$input"
