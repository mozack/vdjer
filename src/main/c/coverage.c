#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <utility>
#include <vector>
#include "quick_map3.h"

using namespace std;

char mate_coverage_is_valid(int read_length, int contig_len, int eval_start, int eval_stop, int read_span,
		               int insert_low, int insert_high, int floor, vector<mapped_pair>& mapped_reads,
		               vector<pair<int,int> >& start_positions, char is_debug) {

//s	floor = 1;
	int num_reads = start_positions.size();
	int begin_idx = 0;
	int pos = eval_start;
	char is_valid = 1;

	int mate_low = 0;
	int mate_high = 0;

	while (mate_low < eval_stop && pos < eval_stop && is_valid) {
		while(begin_idx<num_reads && start_positions[begin_idx].first+read_length-1 < pos) {
			begin_idx += 1;
		}

		mate_low = pos + insert_low - read_length - read_length/2;
		mate_high = pos + insert_high - read_length + read_length/2;

		if (mate_high > eval_stop) {
			mate_high = eval_stop+1;
		}

		int coverage[contig_len+1];
		memset(coverage, 0, sizeof(int)*(contig_len+1));

		int idx = begin_idx;
		while (idx < num_reads && start_positions[idx].first <= pos) {
			for (int i=0; i<read_length; i++) {
				coverage[start_positions[idx].second+i] += 1;
			}
			idx += 1;
		}

		for (int j=mate_low; j<mate_high; j++) {
			if (coverage[j] < floor) {
				is_valid = 0;
				printf("MATE_FLOOR: %s\t%d\t%d\t%d\n", mapped_reads[0].contig_id, pos, j, coverage[j]);
				break;
			}
		}

		pos += 1;
	}

	return is_valid;
}


char coverage_is_valid(int read_length, int contig_len, int eval_start, int eval_stop, int read_span,
		               int insert_low, int insert_high, int floor, vector<mapped_pair>& mapped_reads,
		               vector<pair<int,int> >& start_positions, char is_debug) {


	char is_valid = 1;

	// If read_length=50 and read_span=35, then require floor number of reads every 15 bases.
	int start_gap = read_length - read_span;

	int num_reads = start_positions.size();
	int first_pos_in_range = -1;
	int i = 0;

	while (i < num_reads && start_positions[i].first <= (eval_stop-read_span)+1) {
		if (start_positions[i].first >= eval_start) {

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
				first_pos_in_range = start_positions[i].first;
				if (first_pos_in_range > eval_start+start_gap) {
					// Gap in coverage near eval_start
///					printf("START_GAP @ %d\t%s\n", mapped_reads[i].pos1, mapped_reads[i].contig_id);
					is_valid = 0;
					break;
				}
			}

			// Require at least floor reads starting within start gap of current position
			if (start_positions[i-floor].first < start_positions[i].first-start_gap) {
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
	if (i>floor && start_positions[i-floor].first < start_positions[i].first-start_gap) {
		is_valid = 0;
	}

	// If we've covered the region of interest, coverage should be OK.
	if (i > 0 && num_reads > 0 && start_positions[i].first > eval_stop - read_length && is_valid) {
		is_valid = 1;
	} else {
		if (is_valid = 1) {
			printf("POST_CHECK\n");
		}
//		printf("POST_CHECK. num_reads: %d\t%s\n", num_reads, mapped_reads[i].contig_id);
		is_valid = 0;
	}

	if (is_valid) {
		is_valid = mate_coverage_is_valid(read_length, contig_len, eval_start, eval_stop, read_span,
			   insert_low, insert_high, floor, mapped_reads,
			   start_positions, is_debug);
	}

	return is_valid;
}
