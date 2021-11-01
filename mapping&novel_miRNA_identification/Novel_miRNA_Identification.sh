#!/bin/bash

# mapping
nohup mapper.pl \
/data_3T/yaoxintong/RatBodyMap/raw_fqs/merge320files.fq \
-e -h -i -j \
-k TGGAATTC \
-l 18 -m \
-p /data_3T/yaoxintong/RatBodyMap/reference/UCSC_genome/rn7.fa \
-s ./reads_collapsed.fa \
-t ./reads_vs_refdb.arf \
-v -o 4 &


# novel miRNA identification
nohup miRDeep2.pl \
./reads_collapsed.fa \
/data_3T/yaoxintong/RatBodyMap/reference/UCSC_genome/rn7.fa \
./reads_vs_refdb.arf \
/data_3T/yaoxintong/RatBodyMap/reference/miRBase/rno_mature_ref.fa \
/data_3T/yaoxintong/RatBodyMap/reference/miRBase/mature_other.fa \
/data_3T/yaoxintong/RatBodyMap/reference/miRBase/rno_hairpin_ref.fa \
-t rno 2>report.log &
