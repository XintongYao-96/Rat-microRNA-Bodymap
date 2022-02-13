#!/bin/bash

SampleList_ID=$1
dir_sample_list=/data_3T/yaoxintong/RatBodyMap/scripts_v202201/SampleLists
dir_script=/data_3T/yaoxintong/RatBodyMap/scripts_v202201/tasks
dir_result=/data_3T/yaoxintong/RatBodyMap/raw_result_202201/


for ID in $(cat ${dir_sample_list}/SampleList_${SampleList_ID}.txt);
do 
    echo ${ID};    
    # bash novel_miRNA_mapping.sh ${ID} 
     bash ${dir_script}/ReadCount.sh ${ID} novel_miRNA
    
done

