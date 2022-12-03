#!/bin/bash

#SBATCH -A e30836
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -t 24:00:00
#SBATCH --mem=80G
# load modules you need to use
module load rsem/1.3.0

input="srr-nums-data-2"

while read -r line
do
    cd /projects/e30836/project/group2-f22/expression/$line

    rsem-calculate-expression -p 10 --bam /projects/e30836/project/group2-f22/alignment/$line/Aligned.toTranscriptome.out.bam /projects/e30836/hw1/hg38.index/STAR/gencode.v28 $line
done < "$input"
