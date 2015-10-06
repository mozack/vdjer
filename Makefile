# Make file for ABRA
# libAbra is invoked from the ABRA java code

SRCDIR=src/main/c
SAMTOOLS=samtools-1.2
HTSLIB=samtools-1.2/htslib-1.2.1

samtools:
	$(MAKE) -C $(SAMTOOLS)

vdj:	samtools
	g++ -g -pthread -I$(SRCDIR) -I$(JAVA_HOME)/include -I$(JAVA_HOME)/include/linux -I$(SAMTOOLS) -I$(HTSLIB)  $(SRCDIR)/assembler2_vdj.c $(SRCDIR)/seq_score.c $(SRCDIR)/vj_filter.c $(SRCDIR)/seq_to_kmer.c $(SRCDIR)/hash_utils.c $(SRCDIR)/bam_read.c $(SRCDIR)/quick_map3.c $(SRCDIR)/coverage.c $(SRCDIR)/status.c  $(SAMTOOLS)/libbam.a $(HTSLIB)/libhts.a -lz -lpthread -o vdj

vdjg:
	g++ -g -pthread -I$(SRCDIR) -I$(JAVA_HOME)/include -I$(JAVA_HOME)/include/linux $(SRCDIR)/assembler2_vdj_greedy.c $(SRCDIR)/seq_score.c -o vdjg

vdjg_o:
	g++ -g -pthread -I$(SRCDIR) -I$(JAVA_HOME)/include -I$(JAVA_HOME)/include/linux $(SRCDIR)/assembler2_vdj_greedy_output_graph.c $(SRCDIR)/seq_score.c -o vdjg_o

vdj_d:
	g++ -g -pthread -I$(SRCDIR) -I$(JAVA_HOME)/include -I$(JAVA_HOME)/include/linux $(SRCDIR)/assembler2_vdj_dump.c $(SRCDIR)/seq_score.c -o vdj_d

vdjk11:
	g++ -g -pthread -I$(SRCDIR) -I$(JAVA_HOME)/include -I$(JAVA_HOME)/include/linux $(SRCDIR)/assembler2_vdj_k11.c $(SRCDIR)/seq_score.c -o vdjk11

vdj25:
	g++ -g -I$(SRCDIR) -I$(JAVA_HOME)/include -I$(JAVA_HOME)/include/linux $(SRCDIR)/assembler2.5_vdj.c $(SRCDIR)/seq_score.c -o vdj2.5

seqd:
	g++ -g $(SRCDIR)/seq_dist.c $(SRCDIR)/seq_to_kmer.c -o seqd

vjf:
	g++ -g -I$(SRCDIR) $(SRCDIR)/vj_filter.c $(SRCDIR)/seq_to_kmer.c $(SRCDIR)/hash_utils.c -o vjf

quickmap:
	g++ -g -I$(SRCDIR) $(SRCDIR)/quick_map2.c $(SRCDIR)/hash_utils.c -o quickmap

findcdr3:
	g++ -g -I$(SRCDIR) $(SRCDIR)/cdr3_win.c $(SRCDIR)/vj_filter.c $(SRCDIR)/seq_to_kmer.c -o findcdr3

bamr:
	g++ -g -I$(SRCDIR) -I$(SAMTOOLS) -I$(HTSLIB) $(SRCDIR)/bam_read.c $(SRCDIR)/hash_utils.c $(SRCDIR)/quick_map3.c $(SAMTOOLS)/libbam.a $(HTSLIB)/libhts.a -lz -lpthread -o bamr
