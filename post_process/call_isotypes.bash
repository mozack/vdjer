#!/bin/bash

SEQUENCE_FASTA=$1

cat $SEQUENCE_FASTA | python /datastore/nextgenout4/seqware-analysis/lmose/vdj/scripts/get_cseq.py > vdj_const.fa
/datastore/nextgenout4/seqware-analysis/lmose/vdj/scripts/isotypes_star.bash
samtools sort vdj_const.bam vdj_const.sort
samtools index vdj_const.sort.bam 
samtools view vdj_const.sort.bam chr14 2> st.err | python /datastore/nextgenout4/seqware-analysis/lmose/vdj/scripts/call_isotypes.py /datastore/nextgenout4/seqware-analysis/lmose/vdj/ref/human/hg38/igh_constant_func.bed > isotypes.txt
