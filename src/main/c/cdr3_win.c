#include "seq_dist.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <sparsehash/sparse_hash_map>
#include <sparsehash/sparse_hash_set>
#include <sparsehash/dense_hash_set>
#include "vj_filter.h"

using namespace std;
using google::dense_hash_set;

int main(int argc, char** argv) {
	fprintf(stderr, "Here...\n");

	char* vjf_v_file = argv[1];
	char* vjf_j_file = argv[2];
	char* contig_file = argv[3];
	int vjf_max_dist = atoi(argv[4]);
	int vjf_min_win = atoi(argv[5]);
	int vjf_max_win = atoi(argv[6]);
	char vjf_j_conserved = argv[7][0];
	int vjf_window_span = atoi(argv[8]);
	char vjf_j_extension = atoi(argv[9]);

//	PRINT_CDR3_INDEX = 1;

	fprintf(stderr, "V: %s\nJ: %s\ncontigs: %s\ndist: %d\nmin_win: %d\nmax_win: %d\nj_cons: %c\nwindow_span: %d\nj_extension: %d\n",
			vjf_v_file, vjf_j_file, contig_file, vjf_max_dist, vjf_min_win, vjf_max_win, vjf_j_conserved, vjf_window_span, vjf_j_extension);

	vjf_init(vjf_v_file, vjf_j_file, vjf_max_dist, vjf_min_win, vjf_max_win,
			vjf_j_conserved, vjf_window_span, vjf_j_extension);

	FILE* fp = fopen(contig_file, "r");

	char contig[100000];

	dense_hash_set<const char*, vjf_hash, vjf_eqstr> windows;

	int i = 0;
	while (fgets(contig, 100000, fp) != NULL) {
		if(contig[0] != '>') {
			vjf_search(contig, windows);
		}
	}
}
