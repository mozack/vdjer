STAR=/datastore/nextgenout4/seqware-analysis/lmose/software/star/STAR_2.4.2a/STAR-STAR_2.4.2a/source/STAR

$STAR \
        --runThreadN 8 \
        --genomeDir /datastore/nextgenout4/seqware-analysis/lmose/ref/hg38/star \
        --readFilesIn vdj_const.fa \
        --outSAMunmapped Within \
        --sjdbGTFfile /datastore/nextgenout4/seqware-analysis/lmose/ref/hg38/hg38.knowngene.gtf \
        --sjdbOverhang 100 \
        --outStd SAM | samtools view -1 -bS - > vdj_const.bam
        2> star1.log
