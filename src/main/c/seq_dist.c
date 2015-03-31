#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "seq_dist.h"

using namespace std;

#define MAX_16_BASES 4294967295
#define MAX_DIST 5

int base_val(char base) {
	int val = -1;

	switch(base) {
		case 'A':
			val = 0;
			break;
		case 'T':
			val = 1;
			break;
		case 'C':
			val = 2;
			break;
		case 'G':
			val = 3;
			break;
		default:
			fprintf(stderr, "Error converting base: %c\n", base);
			exit(-1);
	}

	return val;
}

unsigned long seq_to_int(const char* seq) {
	int val = 0;

	for (int i=0; i<SEQ_LEN; i++) {
		val = val << 2;
		val += base_val(seq[i]);
	}

	return val;
}

int edit_dist(unsigned long i1, unsigned long i2) {
	int dist;
	unsigned long val;

	dist = 0;
	val = i1 ^ i2;

	// Count the number of differing 2 bit bases
	while (val != 0) {
		if ((val & 3) != 0) {
			dist++;
		}
		val = val >> 2;
	}

	return dist;
}

int edit_distance(char* seq1, char* seq2) {
	int i1 = seq_to_int(seq1);
	int i2 = seq_to_int(seq2);

	int dist = edit_dist(i1, i2);
	return dist;
}

void get_kmers(char* input, vector<unsigned long>& kmers) {
	FILE* in = fopen(input, "r");

	char kmer[1024];
	while (fgets(kmer, 1024, in) != NULL) {
		unsigned long k = seq_to_int(kmer);
		kmers.push_back(k);
	}

	fclose(in);
}

void process_kmers(char* input, unsigned long start, unsigned long end) {

	vector<unsigned long> kmers;
	get_kmers(input, kmers);

	// for (unsigned long i=0; i<MAX_16_BASES; i++) {
	for (unsigned long i=start; i<=end; i++) {

		int min_dist = SEQ_LEN + 1;

		vector<unsigned long>::const_iterator kmer;
		for (kmer=kmers.begin(); kmer!=kmers.end(); kmer++) {
			int dist = edit_dist(i, *kmer);
			if (dist < min_dist) {
				min_dist = dist;
			}
		}

		if (min_dist <= MAX_DIST) {
			printf("%lu\t%d\n", i, min_dist);
		}
	}
}

int main(int argc, char** argv) {

	if (argc != 4) {
		fprintf(stderr, "Usage: seqd <input> <start> <end>\n");
		exit(-1);
	}

	char* input = argv[1];
	unsigned long start = strtoul(argv[2], NULL, 10);
	unsigned long end   = strtoul(argv[3], NULL, 10);

	process_kmers(input, start, end);

	/*
	char* s1 = "GATTCCGGCGCCGATC";
//	char* s2 = "GATACCGGCGCCGATG";
	char* s2 = "AGAAAAAAAAAAACAA";
	int dist = edit_distance(s1, s2);

	printf("dist: %d\n", dist);

	*/
}
