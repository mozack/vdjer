#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <set>
#include <vector>
#include <sparsehash/dense_hash_map>
#include "seq_score.h"
#include "hash_utils.h"

using namespace std;
using google::dense_hash_map;

char* seq_score_matrix;
int seq_score_size;
char** ref_contigs;
vector<char*> seq_score_contigs;

#define CONTIG_BUF_MAX 10000000
vector<char*> score_seq_contigs;

dense_hash_map<const char*, vector<int>, vregion_hash, vregion_eqstr> contig_index;

//int kmer_size = 35;
//int VREGION_KMER_SIZE = 11;

extern int kmer_size;
extern int VREGION_KMER_SIZE;


void score_seq_init(int max_len1, int max_len2) {
	seq_score_size = (max_len1+1) * (max_len2+1);
	seq_score_matrix = (char*) calloc(seq_score_size, sizeof(int));
}

void add_to_index(char* contig) {
	int stop = strlen(contig)-VREGION_KMER_SIZE;
	for (int i=0; i<stop; i++) {
//		printf("i: %d\n", i);
		if (contig_index.find(contig+i) == contig_index.end()) {
			vector<int> positions;
			positions.push_back(i);
			contig_index[contig+i] = positions;
		} else {
			contig_index[contig+i].push_back(i);
		}
	}
}

void score_seq_init(int max_len1, int max_len2, char* fasta) {
	contig_index.set_empty_key(NULL);
	score_seq_init(max_len1, max_len2);

	fflush(stdout);
	char* contig = (char*) calloc(CONTIG_BUF_MAX, sizeof(char));
	FILE* fp = fopen(fasta, "r");

	while (fgets(contig, 1000000, fp) != NULL) {
		if (contig[0] != '>') {
			// Get rid of newline
			contig[strlen(contig)-1] = '\0';
			score_seq_contigs.push_back(contig);
			add_to_index(contig);
		}

		contig = (char*) calloc(CONTIG_BUF_MAX, sizeof(char));
	}
	fflush(stdout);
	fclose(fp);
}

int seq_score_index(int row, int col, int len1) {
	return col*(len1+1) + row;
}

int score_entry(int row, int col, const char* seq1, const char* seq2, int len1, char* matrix) {
	int val1 = matrix[seq_score_index(row,col-1,len1)] + GAP_PENALTY;
	int val2 = matrix[seq_score_index(row-1,col,len1)] + GAP_PENALTY;
	int val3 = matrix[seq_score_index(row-1,col-1,len1)] + (seq1[row-1] == seq2[col-1] ? MATCH_SCORE : MISMATCH_PENALTY);

	int max = val1 > val2 ? val1 : val2;
	max = max > val3 ? max : val3;

	// Do not allow negative scores
//	max = max < 0 ? 0 : max;

	return max;
}

//
// Identifies max local alignment score for seq1 within seq2
int score_seq(const char* seq1, const char* seq2, int len1, int len2, int threshold) {
	memset(seq_score_matrix, 0, seq_score_size * sizeof(int));

	char* matrix = seq_score_matrix;
	for (int col=1; col<=len2; col++) {
		for (int row=1; row<=len1; row++) {
			matrix[seq_score_index(row, col, len1)] = score_entry(row, col, seq1, seq2, len1, matrix);
//			printf("%d,", matrix[seq_score_index(row, col, len1)]);
		}
//		printf("\n");
	}

	int max = -1;
	int row = len1;
	for (int col=0; col<=len2; col++) {
		for (int row=0; row<=len1; row++) {
			int val = matrix[seq_score_index(row, col, len1)];
			if (val >= threshold) {
				return 1;
			}
		}
	}

	return 0;
}

int score_seq(const char* seq1, int seq_len, int threshold) {

	set<int> to_search;
	int stop = seq_len - VREGION_KMER_SIZE;

	// Identify hits in hash index
	for (int i=0; i<stop; i++) {
		if (contig_index.find(seq1+i) != contig_index.end()) {
			vector<int> positions = contig_index[seq1+i];
			for (vector<int>::iterator it=positions.begin(); it!= positions.end(); ++it) {
				to_search.insert(*it);
			}
		}
	}

//	printf("to_search: %d\t", to_search.size());

	int max = -1;
	for (vector<char*>::const_iterator it=score_seq_contigs.begin(); it != score_seq_contigs.end(); ++it) {
		for (set<int>::iterator idx=to_search.begin(); idx != to_search.end(); ++idx) {
			int start_idx = *idx - kmer_size;
			if (start_idx < 0) {
				start_idx = 0;
			}
			if (start_idx >= strlen(*it)-kmer_size*2) {
				start_idx = strlen(*it)-kmer_size*2 - 1;
			}

			int score = score_seq(seq1, *it + start_idx, seq_len, kmer_size*2, threshold);

			if (score == 1) {
				return 1;
			}
		}
	}

	return 0;
}


int score_seq(const char* seq1, int threshold) {
	return score_seq(seq1, kmer_size, threshold);
}

/*
int main(int argc, char* argv[]) {

//	score_seq_init(35, 500,"/home/lmose/dev/vdj/roots/igh_lv_partial_fix.fa");

	score_seq_init(35, 1000000,"/home/lmose/dev/vdj/roots/v_region_rc.fa");

//	char* root = "CATATCTGTAGACAAGAATAAGAATCAGTTCCTCC";
//	char* root2 = "GTCAGCACGGCATATCTGCAGATCAGCAGCCTATA";
//	char* contig = "ATGGACTGGACCTGGAGGATCCTCTTCTTGGTGGCAGCAGCAACAGGTGCCCACTCCCAGGTGCAGCTGGTGCAATCTGGGTCTGAGTTGAAGAAGCCTGGGGCCTCAGTGAAGGTTTCCTGCAAGGCTTCTGGATACACCTTCACTAGCTATGCTATGAATTGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTGGATGGGATGGATCAACACCAACACTGGGAACCCAACGTATGCCCAGGGCTTCACAGGACGGTTTGTCTTCTCCTTGGACACCTCTGTCAGCACGGCATATCTGCAGATCAGCAGCCTAAAGGCTGAGGACACTGCCGTGTATTACTGTGCGAGAGA";

//	int score = score_seq(root2, contig);

//	char* root2 = "GCAACAGGTGCCCACTGCCAGGTTCAATTGGTGCA";
//
//	int score = score_seq(root2);
//
//	printf("score: %d\n", score);


	FILE* fp = fopen("/home/lmose/dev/vdj/roots/roots2.txt", "r");
//	FILE* fp = fopen("/home/lmose/dev/vdj/roots/t1000.txt", "r");
	int index = 0;
	char root[100];
	memset(root,0,sizeof(char)*100);
	while (fgets(root, 99, fp) != NULL) {
		int last = strlen(root)-1;
		if (last >= 0) {
			root[last] = 0;
		}

		int root_score = score_seq(root, 25);

		printf("SCORE\t%d\t%s\n", root_score, root);

//		printf("HAS_KMER\t%d\t%s\n", has_kmer(root, 35), root);

		memset(root,0,sizeof(char)*100);
	}
	fclose(fp);

}


*/
