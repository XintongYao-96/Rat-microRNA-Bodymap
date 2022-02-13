#!/bin/bash
sample_id=$1
ref_type=transcriptome
dir_result=/data_3T/yaoxintong/RatBodyMap/raw_result_202201/${sample_id}

cat ${dir_result}/${ref_type}/${sample_id}.align2${ref_type}.sam \
| grep -v '@' \
| awk '($2==4) ' \
| awk '{a[$10]++} END{for (i in a){printf("%s\t%d\n", i, a[i])}}' \
| sort -n \
> ${dir_result}/${sample_id}.TranscriptUnmapped.ReadCount
#


# Genome unmapped readcount
sample_id=$1
ref_type=Genome
dir_result=/data_3T/yaoxintong/RatBodyMap/raw_result_202201/${sample_id}

cat ${dir_result}/${ref_type}/${sample_id}.align2${ref_type}.sam \
| grep -v '@' \
| awk '($2==4) ' \
| awk '{a[$10]++} END{for (i in a){printf("%s\t%d\n", i, a[i])}}' \
| sort -n \
> ${dir_result}/${sample_id}.Unmapped.ReadCount
#
