#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "params.h"

void set_chain_info(params* p, char* chain) {
	p->v_region = (char*) calloc(32, sizeof(char));
	p->c_region = (char*) calloc(32, sizeof(char));

	// The coordinates specified here are in the context of hg38
	if (strcmp(chain, "IGH") == 0) {
		strcpy(p->v_region, "chr14:105566277-106879844");
		strcpy(p->c_region, "chr14:105566277-105939754");
		p->j_conserved = 'W';
		p->vj_min_win = 10;
		p->vj_max_win = 90;
	} else if (strcmp(chain, "IGL") == 0) {
		strcpy(p->v_region, "chr22:22026076-22922913");
		strcpy(p->c_region, "chr22:22895375-22922913");
		p->j_conserved = 'F';
		p->vj_min_win = 0;
		p->vj_max_win = 60;
	} else if (strcmp(chain, "IGK") == 0) {
		strcpy(p->v_region, "chr2:89851758-90235368");
		strcpy(p->c_region, "chr2:88857361-88857683");
		p->j_conserved = 'F';
		p->vj_min_win = 0;
		p->vj_max_win = 60;
	} else {
		fprintf(stderr, "Invalid chain specified: %s.  Chain must be one of [IGH,IGL,IGK]\n", chain);
		exit(-1);
	}
}

void set_reference_info(params* p, char* ref_dir) {
	p->v_anchors = (char*) calloc(4096, sizeof(char));
	p->j_anchors = (char*) calloc(4096, sizeof(char));
	p->vdj_fasta = (char*) calloc(4096, sizeof(char));
	p->source_sim_file = (char*) calloc(4096, sizeof(char));

	strcpy(p->v_anchors, ref_dir);
	strcat(p->v_anchors, "/v_index");
	strcpy(p->j_anchors, ref_dir);
	strcat(p->j_anchors, "/j_index");
	strcpy(p->vdj_fasta, ref_dir);
	strcat(p->vdj_fasta, "/ig_vdj.fa");
	strcpy(p->source_sim_file, ref_dir);
	strcat(p->source_sim_file, "/v_region.fa");
}

void set_default_params(params* p) {

	memset(p, 0, sizeof(params));

	p->kmer = 35;
	p->min_node_freq = 3;
	p->min_base_quality = 90;
	p->min_contig_score = -5;
	p->anchor_mismatches = 4;
	p->window_span = 486;
	p->j_extension = 162;
	p->read_filter_floor = 1;
	p->vregion_kmer_size = 15;
	p->min_source_homology_score = 30;
	p->filter_read_span = 35;
	p->filter_mate_span = 48;
	p->eval_start = 50;
	p->eval_stop = 390;
	p->threads = 1;
}
void usage() {
	fprintf(stderr, "vdjer \n");
	fprintf(stderr, "--in <input_bam>\n");
	fprintf(stderr, "--chain <IGH|IGK|IGL>\n");
	fprintf(stderr, "--ref-dir </path/to/vdjer/ref/dir>\n");
	fprintf(stderr, "--mf <min node frequency (default: 3)>\n");
	fprintf(stderr, "--mq <min base quality (default: 90)>\n");
	fprintf(stderr, "--rl <read length>\n");
	fprintf(stderr, "--ms <min contig score (default: -5)\n");
	fprintf(stderr, "--t <threads (default: 1)\n");
	fprintf(stderr, "--am <anchor mismatches (default: 4)\n");
	fprintf(stderr, "--miw <min window length between conserved amino acids>\n");
	fprintf(stderr, "--maw <max window length between conserved amino acids>\n");
	fprintf(stderr, "--jc <conserved J amino acid (W|F)\n");
	fprintf(stderr, "--ws <window span (default: 486)\n");
	fprintf(stderr, "--jext <J extension (default: 162)\n");
	fprintf(stderr, "--vr <V region locus (chr:star-stop)>\n");
	fprintf(stderr, "--cr <C region locus>\n");
	fprintf(stderr, "--ins <expected / median insert length>\n");
	fprintf(stderr, "--rf <read filter floor (default: 1)\n");
	fprintf(stderr, "--k <kmer size (default: 35)>\n");
	fprintf(stderr, "--mrs <min source node homology score (default: 30)\n");
	fprintf(stderr, "--rs <read span distance (default: 35)>\n");
	fprintf(stderr, "--rs <mate span distance (default: 48)>\n");
	fprintf(stderr, "--e0 <start position for contig filtering (default: 50)>\n");
	fprintf(stderr, "--e1 <stop position for contig filtering (default: 390)>\n");
}

void print_params(params* p) {

	fprintf(stderr, "%s\t%s\n", "input", p->input_bam);
	fprintf(stderr, "%s\t%d\n", "min node freq", p->min_node_freq);
	fprintf(stderr, "%s\t%d\n", "min base qual", p->min_base_quality);
	fprintf(stderr, "%s\t%d\n", "max read length", p->read_len);
	fprintf(stderr, "%s\t%d\n", "min contig score (log scaled)", p->min_contig_score);
	fprintf(stderr, "%s\t%d\n", "num threads", p->threads);
	fprintf(stderr, "%s\t%s\n", "v anchor file", p->v_anchors);
	fprintf(stderr, "%s\t%s\n", "j anchor file", p->j_anchors);
	fprintf(stderr, "%s\t%d\n", "max anchor mismatches:", p->anchor_mismatches);
	fprintf(stderr, "%s\t%d\n", "min V/J window", p->vj_min_win);
	fprintf(stderr, "%s\t%d\n", "max V/J window", p->vj_max_win);
	fprintf(stderr, "%s\t%d\n", "conserved J AA", p->j_conserved);
	fprintf(stderr, "%s\t%d\n", "total untrimmed contig length", p->window_span);
	fprintf(stderr, "%s\t%d\n", "extension beyond conserved J AA", p->j_extension);
	// Used to extract 15-mers to compare against mapped reads outside of chain's loci
	fprintf(stderr, "%s\t%s\n", "fasta file containing functional V/D/J sequences", p->vdj_fasta);
	fprintf(stderr, "%s\t%s\n", "variable region locus", p->v_region);
	fprintf(stderr, "%s\t%s\n", "constant region locus", p->c_region);
	fprintf(stderr, "%s\t%d\n", "median insert length", p->insert_len);
	fprintf(stderr, "%s\t%d\n", "read coverage floor for contig filtering", p->read_filter_floor);
	fprintf(stderr, "%s\t%d\n", "kmer", p->kmer);
	// Fasta file of V-Region (in correct orientation) used to identify source nodes in assembly
	fprintf(stderr, "%s\t%s\n", "source node similarity file", p->source_sim_file);
	// Kmer size used to seed root homology search
	fprintf(stderr, "%s\t%d\n", "vregion kmer size", p->vregion_kmer_size);
	fprintf(stderr, "%s\t%d\n", "minimum homology score for source nodes", p->min_source_homology_score);
	fprintf(stderr, "%s\t%d\n", "span for read base filtering", p->filter_read_span);
	fprintf(stderr, "%s\t%d\n", "span for mate base filtering", p->filter_mate_span);
	fprintf(stderr, "%s\t%d\n", "start point for contig filtering", p->eval_start);
	fprintf(stderr, "%s\t%d\n", "stop point for contig filtering", p->eval_stop);
}

char file_exists(char* filename) {
	struct stat   buffer;
	return (stat (filename, &buffer) == 0);
}

void validate_params(params* p) {
	int ok = 1;

	if (p->input_bam == NULL) {
		fprintf(stderr, "Input BAM must be specified\n");
		ok = 0;
	} else if (!file_exists(p->input_bam)) {
		fprintf(stderr, "Could not find input file: %s\n", p->input_bam);
		ok = 0;
	}

	if (p->read_len <= 0) {
		fprintf(stderr, "Read length must be specified and > 0\n");
		ok = 0;
	}

	if (p->v_anchors == NULL) {
		fprintf(stderr, "V anchor file must be specified\n");
		ok = 0;
	} else if (!file_exists(p->v_anchors)) {
		fprintf(stderr, "Could not locate v_index file: %s\n", p->v_anchors);
		ok = 0;
	}

	if (p->j_anchors == NULL) {
		fprintf(stderr, "J anchor file must be specified\n");
		ok = 0;
	} else if (!file_exists(p->j_anchors)) {
		fprintf(stderr, "Could not locate j_index file: %s\n", p->j_anchors);
		ok = 0;
	}

	if (p->vdj_fasta == NULL) {
		fprintf(stderr, "VDJ fasta file must be specified\n");
		ok = 0;
	} else if (!file_exists(p->j_anchors)) {
		fprintf(stderr, "Could not locate vdj_fasta file: %s\n", p->vdj_fasta);
		ok = 0;
	}

	if (p->source_sim_file == NULL) {
		fprintf(stderr, "source_sim_file file must be specified\n");
		ok = 0;
	} else if (!file_exists(p->j_anchors)) {
		fprintf(stderr, "Could not locate source_sim_file file: %s\n", p->source_sim_file);
		ok = 0;
	}

	if ((p->j_conserved != 'W' && p->j_conserved != 'F')) {
		fprintf(stderr, "Conserved J AA must be W or F: %c\n", p->j_conserved);
		ok = 0;
	}

	if (p->v_region == NULL) {
		fprintf(stderr, "V region coordinates must be specified\n");
		ok = 0;
	}

	if (p->c_region == NULL) {
		fprintf(stderr, "C region coordinates must be specified\n");
		ok = 0;
	}

	if (p->insert_len <= 0) {
		fprintf(stderr, "insert_len must be specified and > 0\n");
		ok = 0;
	}

	if (!ok) {
		usage();
		exit(-1);
	}

	print_params(p);
}

char parse_params(int argc, char** argv, params* p) {

	set_default_params(p);

	for (int i=1; i<argc-1; i+=2) {
		char* param = argv[i];
		char* value = argv[i+1];


		if (!strcmp(param, "--help")) {
			usage();
			exit(0);
		} else if (!strcmp(param, "--in")) {
			p->input_bam = value;
		} else if (!strcmp(param, "--chain")) {
			set_chain_info(p, value);
		} else if (!strcmp(param, "--ref-dir")) {
			set_reference_info(p, value);
		} else if (!strcmp(param, "--mf")) {
			p->min_node_freq = atoi(value);
		} else if (!strcmp(param, "--mq")) {
			p->min_base_quality = atoi(value);
		} else if (!strcmp(param, "--rl")) {
			p->read_len = atoi(value);
		} else if (!strcmp(param, "--mcs")) {
			p->min_contig_score = atof(value);
		} else if (!strcmp(param, "--t")) {
			p->threads = atoi(value);
		} else if (!strcmp(param, "--vf")) {
			p->v_anchors = value;
		} else if (!strcmp(param, "--jf")) {
			p->j_anchors = value;
		} else if (!strcmp(param, "--am")) {
			p->anchor_mismatches = atoi(value);
		} else if (!strcmp(param, "--miw")) {
			p->vj_min_win = atoi(value);
		} else if (!strcmp(param, "--maw")) {
			p->vj_max_win = atoi(value);
		} else if (!strcmp(param, "--jc")) {
			p->j_conserved = value[0];
		} else if (!strcmp(param, "--ws")) {
			p->window_span = atoi(value);
		} else if (!strcmp(param, "-jext")) {
			p->j_extension = atoi(value);
		} else if (!strcmp(param, "--vdjf")) {
			p->vdj_fasta = value;
		} else if (!strcmp(param, "--vr")) {
			p->v_region = value;
		} else if (!strcmp(param, "--cr")) {
			p->c_region = value;
		} else if (!strcmp(param, "--ins")) {
			p->insert_len = atoi(value);
		} else if (!strcmp(param, "--rf")) {
			p->read_filter_floor = atoi(value);
		} else if (!strcmp(param, "--k")) {
			p->kmer = atoi(value);
		} else if (!strcmp(param, "--rms")) {
			p->source_sim_file = value;
		} else if (!strcmp(param, "--vk")) {
			p->vregion_kmer_size = atoi(value);
		} else if (!strcmp(param, "--mrs")) {
			p->min_source_homology_score = atoi(value);
		} else if (!strcmp(param, "--rs")) {
			p->filter_read_span = atoi(value);
		} else if (!strcmp(param, "--ms")) {
			p->filter_mate_span = atoi(value);
		} else if (!strcmp(param, "--e0")) {
			p->eval_start = atoi(value);
		} else if (!strcmp(param, "--e1")) {
			p->eval_stop = atoi(value);
		}
	}

	validate_params(p);
}





