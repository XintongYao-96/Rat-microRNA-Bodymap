#!/bin/bash
dir_result=/data_3T/yaoxintong/RatBodyMap/raw_result_202201
#dir_readcount=/data_3T/yaoxintong/RatBodyMap/miRNA_readcount
#dir_readstats=/data_3T/yaoxintong/RatBodyMap/readstats
dir_novel_readcount=/data_3T/yaoxintong/RatBodyMap/novel_readcount

# Collect miRNA readcount
for ID in $(ls ${dir_result});
do 
    echo ${ID}
    #cp ${dir_result}/${ID}/${ID}.miRBase.ReadCount ${dir_readcount}/${ID}.miRBase.ReadCount
    #cp ${dir_result}/${ID}/${ID}.ReadStats ${dir_readstats}/${ID}.ReadStats
    cp ${dir_result}/${ID}/${ID}.novel_miRNA.ReadCount ${dir_novel_readcount}/${ID}.novel_miRNA.ReadCount

done

# Collect readstats
#for ID in $(ls ${dir_result});
#do 
  #  echo ${ID}
  # cp ${dir_result}/${ID}/${ID}.ReadStats ${dir_readstats}/${ID}.ReadStats
    
#done


