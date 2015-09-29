#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include "quick_map3.h"

using namespace std;

char coverage_is_valid(int read_length, int contig_len, int eval_start, int eval_stop, int read_span,
		               int insert_low, int insert_high, int floor, vector<mapped_pair>& mapped_reads) {

	char is_valid = 1;

	// If read_length=50 and read_span=35, then require floor number of reads every 15 bases.
	int start_gap = read_length - read_span;

	int num_reads = mapped_reads.size();
	int first_pos_in_range = -1;
	int i = 0;

	while (i < num_reads && mapped_reads[i].pos1 < (eval_stop-read_span)) {
		if (mapped_reads[i].pos1 >= eval_start) {
			if (i < floor) {
				// Invalid
				is_valid = 0;
				break;
			}

			if (first_pos_in_range < 0) {
				first_pos_in_range = mapped_reads[i].pos1;
				if (first_pos_in_range > eval_start+start_gap) {
					// Gap in coverage near eval_start
					is_valid = 0;
					break;
				}
			}

			// Require at least floor reads starting within start gap of current position
			if (mapped_reads[i-floor].pos1 < mapped_reads[i].pos1-start_gap) {
				// Gap in coverage
				is_valid = 0;
				break;
			}
		}

		i++;
	}

	// If we've covered the region of interest, coverage should be OK.
	if (mapped_reads[i].pos1 > eval_stop - read_length && is_valid) {
		is_valid = 1;
	} else {
		is_valid = 0;
	}

	return is_valid;
}
