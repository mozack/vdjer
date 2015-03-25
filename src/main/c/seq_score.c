#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <list>
#include <set>
#include <sparsehash/sparse_hash_map>
//#include <sparsehash/sparse_hash_set>
#include "seq_score.h"

using namespace std;
using google::sparse_hash_map;

char* seq_score_matrix;
int seq_score_size;
char** ref_contigs;

// KMER size used for placement in VDJ regions
int HASH_KMER = 10;

int sw_count = 0;

/*
struct contig_hash_count {
	char* contig;
	int count;
};
*/

// TODO: Move to separate file
#define BIG_CONSTANT(x) (x##LLU)

uint64_t MurmurHash64A_B ( const void * key, int len, uint64_t seed )
{
  const uint64_t m = BIG_CONSTANT(0xc6a4a7935bd1e995);
  const int r = 47;

  uint64_t h = seed ^ (len * m);

  const uint64_t * data = (const uint64_t *)key;
  const uint64_t * end = data + (len/8);

  while(data != end)
  {
    uint64_t k = *data++;

    k *= m;
    k ^= k >> r;
    k *= m;

    h ^= k;
    h *= m;
  }

  const unsigned char * data2 = (const unsigned char*)data;

  switch(len & 7)
  {
  case 7: h ^= uint64_t(data2[6]) << 48;
  case 6: h ^= uint64_t(data2[5]) << 40;
  case 5: h ^= uint64_t(data2[4]) << 32;
  case 4: h ^= uint64_t(data2[3]) << 24;
  case 3: h ^= uint64_t(data2[2]) << 16;
  case 2: h ^= uint64_t(data2[1]) << 8;
  case 1: h ^= uint64_t(data2[0]);
          h *= m;
  };

  h ^= h >> r;
  h *= m;
  h ^= h >> r;

  return h;
}

struct eqstr_B
{
  bool operator()(const char* s1, const char* s2) const
  {
    return (s1 == s2) || (s1 && s2 && strncmp(s1, s2, HASH_KMER) == 0);
  }
};

struct my_hash_B
{
	uint64_t operator()(const char* kmer) const
	{
		return MurmurHash64A_B(kmer, HASH_KMER, 97);
		//return chunk;
	}
};

sparse_hash_map<const char*, set<char*>*, my_hash_B, eqstr_B>* ref_hash;

// Max size of sequence used for root node comparison
#define MAX_ROOT_SEQUENCE_SIZE 10000

void init_score_seq(int max_len1, int max_len2) {
	seq_score_size = (max_len1+1) * (max_len2+1);
	seq_score_matrix = (char*) calloc(seq_score_size, sizeof(int));
}

int seq_score_index(int row, int col, int len1) {
	return col*(len1+1) + row;
}

//TODO: Incorporate Affine gap penalty?
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
int score_seq(const char* seq1, const char* seq2, int len1, int len2) {
	memset(seq_score_matrix, 0, seq_score_size * sizeof(int));

//	printf("seq_score_size: %d\n", seq_score_size);
//	printf("strlen(seq1): %d\n", strlen(seq1));
//	printf("strlen(seq2): %d\n", strlen(seq2));

	char* matrix = seq_score_matrix;
	for (int col=1; col<=len2; col++) {
		for (int row=1; row<=len1; row++) {
//			printf("seq_score_index(row, col, len1): %d\n", seq_score_index(row, col, len1));
			matrix[seq_score_index(row, col, len1)] = score_entry(row, col, seq1, seq2, len1, matrix);
		}
	}

	//TODO: Assess row locality in memory here for performance
	int max = -1;
	int row = len1;
	for (int col=0; col<=len2; col++) {
		for (int row=0; row<=len1; row++) {
			int val = matrix[seq_score_index(row, col, len1)];
			if (val > max) {
				max = val;
			}
		}
	}

	sw_count += 1;

	return max;
}

int score_seq(const char* seq1, const char* seq2) {
	return score_seq(seq1, seq2, strlen(seq1), strlen(seq2));
}

void init_ref_hash() {
	ref_hash = new sparse_hash_map<const char*, set<char*>*, my_hash_B, eqstr_B>();

	ref_hash->set_deleted_key(NULL);

	int index = 0;
	while (ref_contigs[index] != NULL) {
		for (int i=0; i<strlen(ref_contigs[index])-HASH_KMER; i++) {
			set<char*>* contig_set = (*ref_hash)[ref_contigs[index]+i];
			if (contig_set == NULL) {
				contig_set = new set<char*>();
				(*ref_hash)[ref_contigs[index]+i] = contig_set;
			}

			// Add the current contig to the contig set for this kmer
			contig_set->insert(ref_contigs[index]);
		}
		index += 1;
	}
}

char** init_root_scoring(const char* fasta, int HASH_KMER) {

	init_score_seq(HASH_KMER, MAX_ROOT_SEQUENCE_SIZE);

	// Load each fasta seqeunce
	ref_contigs = (char**) calloc(1000, sizeof(char*));
	FILE* fp = fopen(fasta, "r");
	int index = 0;
	ref_contigs[index] = (char*) calloc(MAX_ROOT_SEQUENCE_SIZE, sizeof(char));
	while (fgets(ref_contigs[index], MAX_ROOT_SEQUENCE_SIZE-1, fp) != NULL && index < 1000) {
		// Replace newline with string terminator
		int last_char = strlen(ref_contigs[index]) - 1;
		if (last_char > 0 && ref_contigs[index][0] != '>') {
			ref_contigs[index][last_char] = '\0';
			index++;
			ref_contigs[index] = (char*) calloc(MAX_ROOT_SEQUENCE_SIZE, sizeof(char));
		}
	}

	if (strlen(ref_contigs[index]) == 0) {
		ref_contigs[index] = NULL;
	}

	init_ref_hash();

	return ref_contigs;
}

#define KMER 25
int score_seq(const char* seq) {

	set<char*> contigs_to_search;

	// KMER = assembly kmer size
	// HASH_KMER = VDJ hash lookup size
	for (int i=0; i<KMER-HASH_KMER; i++) {
//		printf("seq: %s i: %d\n", seq+i, i);
		fflush(stdout);
		set<char*>* hash_contigs = (*ref_hash)[seq+i];
//		printf("hash: %x\n", hash_contigs);
		fflush(stdout);

//		printf("hash: %s\n", seq+i);

		if (hash_contigs != NULL) {
			for (set<char*>::const_iterator iter = hash_contigs->begin(); iter != hash_contigs->end(); ++iter ) {
				char* contig = (*iter);
//				printf("\tcontig: %s\n", contig);
				contigs_to_search.insert(contig);
			}
		}
	}

	int max = -1;
	int index = 0;

//	while (ref_contigs[index] != NULL) {
//		char* contig = ref_contigs[index];
	for (set<char*>::const_iterator iter = contigs_to_search.begin(); iter != contigs_to_search.end(); ++iter ) {
		char* contig = (*iter);

		int score = score_seq(seq, contig, KMER, strlen(contig));
//		printf("score: %d\n", score);
		if (score > max) {
			max = score;
		}
		index += 1;
	}

	return max;
}

/*
int main(int argc, char* argv[]) {

//	char** ref_contigs = init_root_scoring("/home/lmose/code/abra/mouse.bcr.fa", 25);
	char** ref_contigs = init_root_scoring("/home/lmose/code/abra/mm10.bcr.constant.fa", 25);

	int idx = 0;
//	while (ref_contigs[idx] != NULL) {
//		printf("contig: [%s]\n", ref_contigs[idx++]);
//	}

////	int score = score_seq("ATTTTATTCGAAAAATCGACGGAGT");
//	int score = score_seq("AAGGCGCCTGTGTCTGGCAGAGCCG");
////	int score = score_seq("CTGGTCCCTGGGGCACGCCTCGGCC");
//
//	printf("Final score: %d\n", score);
//	printf("Num SW searches: %d\n", sw_count);


	FILE* fp = fopen("roots.txt", "r");
	int index = 0;
	char root[100];
	memset(root,0,sizeof(char)*100);
	while (fgets(root, 99, fp) != NULL) {
		int last = strlen(root)-1;
		if (last >= 0) {
			root[last] = 0;
		}

		int root_score = score_seq(root);

		printf("SCORE\t%d\t%s\n", root_score, root);
		memset(root,0,sizeof(char)*100);
	}
	fclose(fp);


//	init_score_seq(100, 100);
//
//	int score = score_seq("ATCGA", "ATTTTATTCGAAAAA");
//
//	printf("score: %d\n", score);
//
//	score = score_seq("TTTTTTT", "ATTTTATTCGAAAAA");
//
//	printf("score: %d\n", score);
//
//	score = score_seq("CCGGCCA", "ATTTTATTCGAAAAA");
//
//	printf("score: %d\n", score);
}
*/
