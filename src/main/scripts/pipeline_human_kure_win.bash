#!/bin/bash

export PATH=/datastore/nextgenout4/seqware-analysis/lmose/software/bwa/bwa-0.7.9a:/datastore/nextgenout4/seqware-analysis/lmose/software/samtools-0.1.19:$PATH

TYPE=$1
UUID=$2

cd /datastore/nextgenout4/seqware-analysis/lmose/vdj/breast/$TYPE/$UUID

rm mapsplice.sort.bam
rm mapsplice.sort.bam.bai

samtools sort -@ 8 -m 2G mapsplice_out/alignments.bam mapsplice.sort

samtools index mapsplice.sort.bam

#rm -rf vdj_analysis2.bak

mv vdj_analysis2 vdj_analysis2.bak
mkdir vdj_analysis2
cd vdj_analysis2

# Extract reads from VDJ regions, C regions, their mates and unmapped read pairs where either read in the pair contains a 15-mer in VDJ regions.  Reads are extracted into a BAM file
java -Xmx8G -cp /datastore/nextgenout4/seqware-analysis/lmose/vdj/mouse/abra/target/abra-0.87-SNAPSHOT-jar-with-dependencies.jar vdj.VdjReadExtractor2 ../mapsplice.sort.bam vdj.reads.bam /datastore/nextgenout4/seqware-analysis/lmose/vdj/human/igh_vdj_no_p.bed /datastore/nextgenout4/seqware-analysis/lmose/vdj/human/igh_constant_no_p.bed /datastore/nextgenout4/seqware-analysis/lmose/vdj/human/vdj_no_p.fa 15 > extract.log 2>&1

# Convert BAM file into vdj assembler input format

java -Xmx2G -cp /datastore/nextgenout4/seqware-analysis/lmose/vdj/mouse/abra/target/abra-0.87-SNAPSHOT-jar-with-dependencies.jar vdj.ReadExtractorFull2 vdj.reads.bam vdj.reads > extract2.log 2>&1

# Assemble - produces vdj.fa  - Currently hardcoded to use 8 threads in code
/datastore/nextgenout4/seqware-analysis/lmose/vdj/mouse/samples/vdj.bash vdj.reads 10

#mv get_cdr3.log get_cdr3.log.bak
#rm vdj.analysis.finis
#rm cdr3.txt cdr3.fa cdr3.reads.bam cdr3.scores.txt cdr3.best.txt vdj.prediction vdj.prediction.aa

# Extract CDR3 windows
java -Xmx8G -cp /datastore/nextgenout4/seqware-analysis/lmose/vdj/mouse/abra/target/abra-0.87-SNAPSHOT-jar-with-dependencies.jar vdj.GetCdr3Windows2 vdj.fa cdr3.txt 90 > get_cdr3.log 2>&1 

# Convert cdr file to fasta format
cat cdr3.txt | cut -f 2- | tr '\t' '\n' | cat -n | awk '{print ">cdr3_" $1 "\n" $2}' > cdr3.fa

# Map original reads to CDR3 windows
bwa index cdr3.fa
#cat ../1.fastq ../2.fastq > both.fastq
#cat ../*.fastq > both.fastq
cat ../*.fastq | bwa mem -t 8 cdr3.fa - 2> bwa.log | samtools view -F 0x04 -bS - > cdr3.reads.bam

# CDR3 HMM scoring
java -Xmx16G -cp /datastore/nextgenout4/seqware-analysis/lmose/vdj/mouse/hmm/abra-0.87-SNAPSHOT-jar-with-dependencies.jar:/datastore/nextgenout4/seqware-analysis/lmose/vdj/mouse/hmm/GenomeAnalysisTK.jar vdj.hmm.ContigSelector cdr3.reads.bam cdr3.txt cdr3.scores.txt 8 > cdr3.score.log

# ID top scoring CDR3
cat cdr3.scores.txt | sort -n -k1,1 | tail -1 > cdr3.best.txt

# Extend CDR3 window on 5' and 3' ends
cat cdr3.best.txt | /datastore/nextgenout4/seqware-analysis/lmose/software/python2.7.8/Python-2.7.8/python /datastore/nextgenout4/seqware-analysis/lmose/vdj/mouse/samples/extend_cdr3.py > vdj.prediction 2> vdj.prediction.log

# Extract unique AA from vdj.prediction(s)
cut -f 1 vdj.prediction | sort | uniq > vdj.prediction.aa

touch vdj.analysis.finis
