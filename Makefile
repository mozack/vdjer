# Make file for V'DJer

SRCDIR=src/main/c
SAMTOOLS=samtools-1.2
HTSLIB=samtools-1.2/htslib-1.2.1

samtools:
	$(MAKE) -C $(SAMTOOLS)

vdj:	samtools
	g++ -g -pthread -I$(SRCDIR) -I$(JAVA_HOME)/include -I$(JAVA_HOME)/include/linux -I$(SAMTOOLS) -I$(HTSLIB)  $(SRCDIR)/assembler2_vdj.c $(SRCDIR)/seq_score.c $(SRCDIR)/vj_filter.c $(SRCDIR)/seq_to_kmer.c $(SRCDIR)/hash_utils.c $(SRCDIR)/bam_read.c $(SRCDIR)/quick_map3.c $(SRCDIR)/coverage.c $(SRCDIR)/status.c  $(SAMTOOLS)/libbam.a $(HTSLIB)/libhts.a -lz -lpthread -o vdj

seqd:
	g++ -g $(SRCDIR)/seq_dist.c $(SRCDIR)/seq_to_kmer.c -o seqd

quickmap:
	g++ -g -I$(SRCDIR) $(SRCDIR)/quick_map2.c $(SRCDIR)/hash_utils.c -o quickmap

bamr:
	g++ -g -I$(SRCDIR) -I$(SAMTOOLS) -I$(HTSLIB) $(SRCDIR)/bam_read.c $(SRCDIR)/hash_utils.c $(SRCDIR)/quick_map3.c $(SAMTOOLS)/libbam.a $(HTSLIB)/libhts.a -lz -lpthread -o bamr

ss:
	g++ -g -I$(SRCDIR) $(SRCDIR)/seq_score.c $(SRCDIR)/hash_utils.c -o ss
