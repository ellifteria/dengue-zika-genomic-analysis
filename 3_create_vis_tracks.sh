#!/bin/bash

#SBATCH -A e30836
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -t 24:00:00
#SBATCH --mem=80G

module load samtools/1.6
module load deeptools/3.1.1

input="srr-nums"

while read -r line
do
	mkdir /projects/e30836/project/group2-f22/alignment/$line/track
	cd /projects/e30836/project/group2-f22/alignment/$line/track
	
	samtools view -b ../Aligned.out.sam > ../Aligned.out.bam
	samtools sort ../Aligned.out.bam > Aligned.out.sort.bam
	samtools index Aligned.out.sort.bam
	bamCoverage --bam Aligned.out.sort.bam --normalizeUsing CPM --outFileName "$line".-.bigWig --filterRNAstrand reverse --binSize 1 --numberOfProcessors 60
	bamCoverage --bam Aligned.out.sort.bam --normalizeUsing CPM --outFileName "$line".+.bigWig --filterRNAstrand forward --binSize 1 --numberOfProcessors 60
done < "$input"
