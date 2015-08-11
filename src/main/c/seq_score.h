#ifndef __SEQ_SCORE__
#define __SEQ_SCORE__ 1

#define GAP_PENALTY -2
#define MATCH_SCORE 1
#define MISMATCH_PENALTY -1

void init_score_seq(int max_len1, int max_len2);

int score_seq(const char* seq1, const char* seq2, int len1, int len2);

char** init_root_scoring(const char* fasta, int kmer_size);

// Score the input sequence against all sequences loaded via init_root_scoring
// returning the max score.
int score_seq(const char* seq);

#endif
