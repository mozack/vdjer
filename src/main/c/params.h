#ifndef __PARAMS__
#define __PARAMS__

struct params {
	char* input_bam;
	int min_node_freq;
	int min_base_quality;
	float min_contig_score;
	int threads;
	char* v_anchors;
	char* j_anchors;
	int anchor_mismatches;
	int vj_min_win;
	int vj_max_win;
	int j_conserved;
	int window_span;
	int j_extension;
	char* vdj_fasta;
	char* v_region;
	char* c_region;
	int insert_len;
	int read_filter_floor;
	int kmer;
	char* source_sim_file;
	int vregion_kmer_size;
	int min_source_homology_score;
	int filter_read_span;
	int filter_mate_span;
	int eval_start;
	int eval_stop;
	int window_overlap_check_size;
};

char parse_params(int argc, char** argv, params* p);

#endif
