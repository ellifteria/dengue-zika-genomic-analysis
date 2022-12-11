#!/bin/bash

#SBATCH -A e30836
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -t 24:00:00
#SBATCH --mem=80G

module load sratoolkit/3.0.0

input="./srr-nums"

while read -r line
do
	fastq-dump -I --split-files $line -O /projects/e30836/project/group2-f22/fastq-data/$line
done < "$input"
