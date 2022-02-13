#!/bin/bash

sample_id=$1

# Project directory
dir_script=/data_3T/yaoxintong/RatBodyMap/scripts_v202201/tasks
# mkdir /data_3T/yaoxintong/RatBodyMap/result_202201/${sample_id}
dir_result=/data_3T/yaoxintong/RatBodyMap/raw_result_202201/${sample_id}

### Map to novel miRNA
mkdir ${dir_result}/novel_miRNA
in_fastq=${dir_result}/TrimAdapter/${sample_id}.trimAdapt.filter.fastq
dir_output=${dir_result}/novel_miRNA
refname=novel_miRNA
dir_index=/data_3T/yaoxintong/RatBodyMap/reference/novel_miRNA/rno_mature_novel.fa

bash ${dir_script}/Mapping_RNA.sh \
    ${sample_id} \
    ${in_fastq} \
    ${dir_output} \
    ${refname} \
    ${dir_index} 