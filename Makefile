# Make file for ABRA
# libAbra is invoked from the ABRA java code

SRCDIR=src/main/c

vdj:
	g++ -g -pthread -I$(SRCDIR) -I$(JAVA_HOME)/include -I$(JAVA_HOME)/include/linux $(SRCDIR)/assembler2_vdj.c $(SRCDIR)/seq_score.c -o vdj

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
	g++ -g $(SRCDIR)/seq_dist.c -o seqd

vjf:
	g++ -g -I$(SRCDIR) $(SRCDIR)/vj_filter.c -o vjf
