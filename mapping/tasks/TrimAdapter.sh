#!/bin/bash

### input
sample_id=$1
dir_result=$2

### Project directory
# output directory
cd ${dir_result}
mkdir ./TrimAdapter
cd ./TrimAdapter

#input directory
dir_unzip_data=/data_3T/yaoxintong/RatBodyMap/data/unzip_raw_data
in_fastq=${dir_unzip_data}/${sample_id}.fq


## Clipping adapter
adapter_seq=TGGAATTCTCGGGTGCCAAGG

fastp -Q \
    --length_required 0\
    --adapter_sequence ${adapter_seq} \
    -i ${in_fastq} \
    -o ./${sample_id}.trimAdapt.fastq.gz \
    2>> ./${sample_id}.trimAdapt.log

## Length Distribution
echo -e "Length\tReadCount" > $1.trimAdapt.lengthDistribute
      zcat $1.trimAdapt.fastq.gz | paste - - - - | cut -f 2 | 
            awk '{a[length($1)]++}END{for(i in a){print i,a[i]}}' | sort -n \
            >> ${sample_id}.trimAdapt.lengthDistribute


## Filtering
qualified_quality_phred=20
unqualified_percent_limit=20
n_base_limit=2
length_required=16
length_limit=35

fastp -A \
      --length_required ${length_required} \
      --length_limit ${length_limit} \
      --qualified_quality_phred ${qualified_quality_phred} \
      --unqualified_percent_limit ${unqualified_percent_limit} \
      --n_base_limit ${n_base_limit} \
      -i ${sample_id}.trimAdapt.fastq.gz \
      -o ${sample_id}.trimAdapt.filter.fastq.gz \
      2> ${sample_id}.filter.log







