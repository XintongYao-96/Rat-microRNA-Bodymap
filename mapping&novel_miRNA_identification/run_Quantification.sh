#!/bin/bash
for ID in $(cat SampleList_$1.txt);
do 
    echo ${ID};    
    bash Quantification.sh ${ID} > ${dir_result_log}/${ID}.log;
done
