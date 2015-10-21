#ifndef __SEQ_SCORE__
#define __SEQ_SCORE__ 1

// Sequence similarity scoring
#define GAP_PENALTY -1
#define MATCH_SCORE 1
//#define MISMATCH_PENALTY -1
#define MISMATCH_PENALTY 0

//
// Initialize scoring matrix and load contigs from input fasta
void score_seq_init(int max_len1, int max_len2, char* fasta);

// Score the input sequence against all sequences loaded via init_root_scoring
// returning the max score.
int score_seq(const char* seq, int seq_len);

#endif
