#ifndef __VJ_FILTER__
#define __VJ_FILTER__

//char PRINT_CDR3_INDEX = 0;

extern __thread char* vjf_cdr3_block_buffer;

// Init V and J anchor indices as well as search params
void vjf_init(char* v_file, char* j_file, int max_dist, int min_win, int max_win,
		char j_conserved, int window_span, int j_extension);

// Search for candidate VDJ windows
void vjf_search(char* contig, google::dense_hash_set<const char*, vjf_hash, vjf_eqstr>& windows);

#endif
