#!/bin/bash

SampleList_ID=$1
dir_sample_list=/data_3T/yaoxintong/RatBodyMap/scripts_v202201/SampleLists
dir_scripts=/data_3T/yaoxintong/RatBodyMap/scripts_v202201

for ID in $(cat ${dir_sample_list}/SampleList_${SampleList_ID}.txt);
do 
    echo ${ID};    
    # bash Workflow.sh ${ID} 
    bash ./RNA_unmapped_count.sh ${ID}
done




