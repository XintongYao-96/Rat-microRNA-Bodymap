#!/bin/bash

sample_id=$1

# Project directory
dir_script=/data_3T/yaoxintong/RatBodyMap/scripts_v202201/tasks
mkdir /data_3T/yaoxintong/RatBodyMap/result_202201/${sample_id}
dir_result=/data_3T/yaoxintong/RatBodyMap/result_202201/${sample_id}

# --- Begin!-----------------------------------------------------------
### Trim Adapter
bash ${dir_script}/TrimAdapter.sh ${sample_id} ${dir_result}


### Map to miRBase
mkdir ${dir_result}/miRBase
in_fastq=${dir_result}/TrimAdapter/${sample_id}.trimAdapt.filter.fastq
dir_output=${dir_result}/miRBase
refname=miRBase
dir_index=/data_3T/yaoxintong/RatBodyMap/reference/miRBase_expand/mature_expand/rno_miRBase_expand.fa

bash ${dir_script}/Mapping_RNA.sh \
    ${sample_id} \
    ${in_fastq} \
    ${dir_output} \
    ${refname} \
    ${dir_index} 


### Map to piRBase
mkdir ${dir_result}/piRBase
in_fastq=${dir_result}/miRBase/${sample_id}.miRBaseUnaligned.fastq
dir_output=${dir_result}/piRBase
refname=piRBase
dir_index=/data_3T/yaoxintong/RatBodyMap/reference/piRBase/rno_piRBase.fa

bash ${dir_script}/Mapping_RNA.sh \
    ${sample_id} \
    ${in_fastq} \
    ${dir_output} \
    ${refname} \
    ${dir_index} 


### Map to GtRNA
mkdir ${dir_result}/GtRNAdb
in_fastq=${dir_result}/piRBase/${sample_id}.piRBaseUnaligned.fastq
dir_output=${dir_result}/GtRNAdb
refname=GtRNAdb
dir_index=/data_3T/yaoxintong/RatBodyMap/reference/GtRNAdb/rn7-tRNAs.fa

bash ${dir_script}/Mapping_RNA.sh \
    ${sample_id} \
    ${in_fastq} \
    ${dir_output} \
    ${refname} \
    ${dir_index} 

### Map to transcriptome
mkdir ${dir_result}/transcriptome
in_fastq=${dir_result}/GtRNAdb/${sample_id}.GtRNAdbUnaligned.fastq
dir_output=${dir_result}/transcriptome
refname=transcriptome
dir_index=/data_3T/yaoxintong/RatBodyMap/reference/NCBI_transcriptome/rno_transcriptome_test/rno_mRNA_test.fa

bash ${dir_script}/Mapping_RNA.sh \
    ${sample_id} \
    ${in_fastq} \
    ${dir_output} \
    ${refname} \
    ${dir_index} 


### Map the Unaligned reads to Rat genome
mkdir ${dir_result}/Genome
in_fastq=${dir_result}/transcriptome/${sample_id}.transcriptomeUnaligned.fastq
dir_output=${dir_result}/Genome
refname=Genome
dir_index=/data_3T/yaoxintong/RatBodyMap/reference/UCSC_genome/rn7.fa

bash ${dir_script}/Mapping_Genome.sh \
    ${sample_id} \
    ${in_fastq} \
    ${dir_output} \
    ${refname} \
    ${dir_index} \


### Map the total filtered reads to Rat Genome
mkdir ${dir_result}/Genome_total
in_fastq=${dir_result}/TrimAdapter/${sample_id}.trimAdapt.filter.fastq
dir_output=${dir_result}/Genome_total
refname=Genome
dir_index=/data_3T/yaoxintong/RatBodyMap/reference/UCSC_genome/rn7.fa

bash ${dir_script}/Mapping_Genome.sh \
    ${sample_id} \
    ${in_fastq} \
    ${dir_output} \
    ${refname} \
    ${dir_index} \



### Reads Statistics
bash ${dir_script}/ReadStats.sh ${sample_id} 



### miRNA expressions
 bash ${dir_script}/ReadCount.sh ${sample_id} miRBase
 bash ${dir_script}/ReadCount.sh ${sample_id} piRBase
 bash ${dir_script}/ReadCount.sh ${sample_id} GtRNAdb
 bash ${dir_script}/ReadCount.sh ${sample_id} transcriptome