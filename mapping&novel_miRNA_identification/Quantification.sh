#!/bin/bash

## Project directory
dir_zip_data=/data_3T/yaoxintong/RatBodyMap/data/zip_raw_data
dir_unzip_data=/data_3T/yaoxintong/RatBodyMap/data/unzip_raw_data
dir_result_mapper=/data_3T/yaoxintong/RatBodyMap/result_mapping_320samples/mapper
dir_result_quantification=/data_3T/yaoxintong/RatBodyMap/result_mapping_320samples/quantification

## Gunzip fq.gz
gunzip -c ${dir_zip_data}/$1.fq.gz > ${dir_unzip_data}/$1.fq

## Map to Genome
mkdir ${dir_result_mapper}/$1
cd ${dir_result_mapper}/$1

mapper.pl \
${dir_unzip_data}/$1.fq \
-e -h -i -j \
-k TGGAATTC \
-l 18 -m \
-p /data_3T/yaoxintong/RatBodyMap/reference/UCSC_genome/rn7.fa \
-s ./$1_reads_collapsed.fa \
-t ./$1_reads_vs_refdb.arf \
-v -o 4


# quantification
mkdir ${dir_result_quantification}/$1
cd ${dir_result_quantification}/$1

quantifier.pl \
-p /data_3T/yaoxintong/RatBodyMap/reference/merge_mature_hairpin/rno_hairpin_merge.fa \
-m /data_3T/yaoxintong/RatBodyMap/reference/merge_mature_hairpin/rno_mature_merge.fa \
-r ${dir_result_mapper}/$1/$1_reads_collapsed.fa \
-t rno \
-y $1
