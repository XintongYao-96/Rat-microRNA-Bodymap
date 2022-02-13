#!/bin/bash
# input
sample_id=$1
in_fastq=$2
dir_output=$3
refname=$4
dir_index=$5
max_mismatch_allowed=2

# unzip
gunzip ${in_fastq}.gz

# align
 bowtie --threads 6 \
       ${dir_index} \
       -v ${max_mismatch_allowed} \
       -q ${in_fastq} \
       -S ${dir_output}/${sample_id}.align2${refname}.sam \
       2> ${dir_output}/${sample_id}.align2${refname}.log   

gzip ${in_fastq}