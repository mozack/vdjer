#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include "quick_map3.h"

using namespace std;

char coverage_is_valid(int read_length, int contig_len, int eval_start, int eval_stop, int read_span,
		               int insert_low, int insert_high, int floor, vector<mapped_pair>& mapped_reads,
		               vector<int>& start_positions, char is_debug) {


	char is_valid = 1;

	// If read_length=50 and read_span=35, then require floor number of reads every 15 bases.
	int start_gap = read_length - read_span;

	int num_reads = start_positions.size();
	int first_pos_in_range = -1;
	int i = 0;

	while (i < num_reads && start_positions[i] < (eval_stop-read_span)) {
		if (start_positions[i] >= eval_start) {

//			if (is_debug) {
//				printf("%s\t%d\t%d\n", mapped_reads[i].contig_id, i, mapped_reads[i].pos1);
//			}

			if (i < floor) {
				// Invalid
//				printf("BEGIN_FLOOR @ %d\t%s\n", mapped_reads[i].pos1, mapped_reads[i].contig_id);
				is_valid = 0;
				break;
			}

			if (first_pos_in_range < 0) {
				first_pos_in_range = start_positions[i];
				if (first_pos_in_range > eval_start+start_gap) {
					// Gap in coverage near eval_start
///					printf("START_GAP @ %d\t%s\n", mapped_reads[i].pos1, mapped_reads[i].contig_id);
					is_valid = 0;
					break;
				}
			}

			// Require at least floor reads starting within start gap of current position
			if (start_positions[i-floor] < start_positions[i]-start_gap) {
				// Gap in coverage
//				printf("MID_FLOOR @ %d\t%s\n", mapped_reads[i].pos1, mapped_reads[i].contig_id);
				is_valid = 0;
				break;
			}
		}

		i++;
	}

	if (i >= num_reads) {
		i = num_reads - 1;
	}

	// Check for gap at last processed read
	if (i>floor && start_positions[i-floor] < start_positions[i]-start_gap) {
		is_valid = 0;
	}

	// If we've covered the region of interest, coverage should be OK.
	if (i > 0 && num_reads > 0 && start_positions[i] > eval_stop - read_length && is_valid) {
		is_valid = 1;
	} else {
		if (is_valid = 1) {
			printf("POST_CHECK\n");
		}
//		printf("POST_CHECK. num_reads: %d\t%s\n", num_reads, mapped_reads[i].contig_id);
		is_valid = 0;
	}

	return is_valid;
}
