#!/bin/bash

# Transcriptome unmapped readcount
sample_id=$1
ref_type=transcriptome
dir_result=/data_3T/yaoxintong/RatBodyMap/raw_result_202201/${sample_id}
dir_collect_result=/data_3T/yaoxintong/RatBodyMap/RNAunmapped_readcount

cat ${dir_result}/${ref_type}/${sample_id}.align2${ref_type}.sam \
| grep -v '@' \
| awk '($2==4) ' \
| awk '{a[$10]++} END{for (i in a){printf("%s\t%d\n", i, a[i])}}' \
| sort -n \
> ${dir_result}/${sample_id}.TranscriptUnmapped.ReadCount

cp ${dir_result}/${sample_id}.TranscriptUnmapped.ReadCount ${dir_collect_result}/${sample_id}.TranscriptUnmapped.ReadCount



