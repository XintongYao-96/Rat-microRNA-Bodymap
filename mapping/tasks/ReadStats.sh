#!/bin/bash
sample_id=$1
dir_input=/data_3T/yaoxintong/RatBodyMap/result_202201/${sample_id}


###### Count reads of each class
n_total_sequence=$(cat ${dir_input}/TrimAdapter/${sample_id}.trimAdapt.log| grep 'total reads' | head -n 1 | cut -d ':' -f 2 | sed 's/ //g')
n_forMapping=$(cat ${dir_input}/TrimAdapter/${sample_id}.filter.log| grep 'reads passed filter' | head -n 1 | cut -d ':' -f 2 | sed 's/ //g')
n_miRNA=$(cat ${dir_input}/miRBase/${sample_id}.align2miRBase.log| grep 'reads with at least one alignment' | head -n 1 | cut -d ':' -f 2 | sed 's/ //g' | cut -d '(' -f 1)
n_piRNA=$(cat ${dir_input}/piRBase/${sample_id}.align2piRBase.log| grep 'reads with at least one alignment' | head -n 1 | cut -d ':' -f 2 | sed 's/ //g' | cut -d '(' -f 1)
n_tRNA=$(cat ${dir_input}/GtRNAdb/${sample_id}.align2GtRNAdb.log| grep 'reads with at least one alignment' | head -n 1 | cut -d ':' -f 2 | sed 's/ //g' | cut -d '(' -f 1)
n_otherGenome=$(cat ${dir_input}/Genome/${sample_id}.align2Genome.log| grep 'reads with at least one alignment' | head -n 1 | cut -d ':' -f 2 | sed 's/ //g' | cut -d '(' -f 1)
n_unmapped=$(cat ${dir_input}/Genome/${sample_id}.align2Genome.log| grep 'reads that failed to align' | head -n 1 | cut -d ':' -f 2 | sed 's/ //g' | cut -d '(' -f 1)
n_totalGenome=$(cat ${dir_input}/Genome_total/${sample_id}.align2Genome.log| grep 'reads with at least one alignment' | head -n 1 | cut -d ':' -f 2 | sed 's/ //g' | cut -d '(' -f 1)


###### Transcriptome subclass
## Assigned reads mapped to transcriptome to each RNA subclass
cat ${dir_input}/transcriptome/${sample_id}.align2transcriptome.sam \
    | grep -v '@' \
    | awk '($2!=4)' \
    | cut -f 3 \
    | sed 's/,_/\t/g' \
    | awk '{print $NF}' \
    > ${dir_input}/transcriptome/${sample_id}.RNA.subcalss.txt

## Count
n_transcriptome=$(cat ${dir_input}/transcriptome/${sample_id}.align2transcriptome.log| grep 'reads with at least one alignment' | head -n 1 | cut -d ':' -f 2 | sed 's/ //g' | cut -d '(' -f 1)
n_mRNA=$(cat ${dir_input}/transcriptome/${sample_id}.RNA.subcalss.txt | grep 'mRNA' | wc -l) 
n_rRNA=$(cat ${dir_input}/transcriptome/${sample_id}.RNA.subcalss.txt | grep -E 'ribosomal_RNA|rRNA' | wc -l)
n_misc_RNA=$(cat ${dir_input}/transcriptome/${sample_id}.RNA.subcalss.txt | grep 'misc_RNA' | wc -l)
n_other_RNA=$((${n_transcriptome}-${n_mRNA}-${n_rRNA}-${n_misc_RNA}))




###### Write ta
file_output=${dir_input}/${sample_id}.ReadStats
echo -e "Type\tReadCount" >> $file_output
echo -e "Total_sequence\t${n_total_sequence}" >> $file_output
echo -e "Total_forAlign\t${n_forMapping}" >> $file_output
echo -e "total_genome\t${n_totalGenome}" >> $file_output
echo -e "miRNA\t${n_miRNA}" >> $file_output
echo -e "piRNA\t${n_piRNA}" >> $file_output
echo -e "mRNA\t${n_mRNA}" >> $file_output
echo -e "tRNA\t${n_tRNA}" >> $file_output
echo -e "rRNA\t${n_rRNA}" >> $file_output
echo -e "miscRNA\t${n_misc_RNA}" >> $file_output
echo -e "other_RNA\t${n_other_RNA}" >> $file_output
echo -e "other_from_Genome\t${n_otherGenome}" >> $file_output
echo -e "unmapped\t${n_unmapped}" >> $file_output


