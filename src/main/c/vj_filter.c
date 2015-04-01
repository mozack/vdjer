#include "seq_dist.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <sparsehash/sparse_hash_map>
#include <sparsehash/sparse_hash_set>

using namespace std;
using google::sparse_hash_set;

// Rough defaults with vmer and jmer 10 bases from edges of annotated v and j genes
#define MIN_WINDOW 10
#define MAX_WINDOW 110

struct eqkmer
{
  bool operator()(unsigned long l1, unsigned long l2) const
  {
    return l1 == l2;
  }
};

#define BIG_CONSTANT(x) (x##LLU)
uint64_t MurmurHash64A ( const void * key, int len, uint64_t seed )
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

struct kmer_hash
{
	uint64_t operator()(unsigned long kmer) const
	{
		return MurmurHash64A(&kmer, 4, 97);
		//return chunk;
	}
};

sparse_hash_set<unsigned long, kmer_hash, eqkmer>* vmers;
sparse_hash_set<unsigned long, kmer_hash, eqkmer>* jmers;


void load_kmers(char* input, sparse_hash_set<unsigned long, kmer_hash, eqkmer>* kmers, int max_dist) {
	FILE* in = fopen(input, "r");

	unsigned long kmer;
	int freq;
	while (fscanf(in, "%lu\t%d\n", &kmer, &freq) == 2) {
//		printf("read: %lu, %d\n", kmer, freq);
		if (freq <= max_dist) {
			kmers->insert(kmer);
		}
	}

	fclose(in);
}

char matches_vmer(unsigned long kmer) {
	sparse_hash_set<unsigned long, kmer_hash, eqkmer>::const_iterator it = vmers->find(kmer);
	return it != vmers->end();
}

char matches_jmer(unsigned long kmer) {
	sparse_hash_set<unsigned long, kmer_hash, eqkmer>::const_iterator it = jmers->find(kmer);
	return it != jmers->end();
}

void print_windows(char* contig) {
	char* contig_index = contig;
	vector<int> v_indices;
	vector<int> j_indices;

	int len = strlen(contig) - SEQ_LEN;
	for (int i=0; i<len; i++) {
		unsigned long kmer = seq_to_int(contig_index);
		if (matches_vmer(kmer)) {
			v_indices.push_back(i);
		}

		if (matches_jmer(kmer)) {
			j_indices.push_back(i);
		}

		contig_index += 1;
	}

	// TODO: traverse vectors in parallel and more intelligently.
	//       no need to compare all values
	vector<int>::const_iterator v;
	for (v=v_indices.begin(); v!=v_indices.end(); v++) {
		vector<int>::const_iterator j;
		for (j=j_indices.begin(); j!=j_indices.end(); j++) {
			int window = *j - *v - SEQ_LEN;

			if (window >= MIN_WINDOW && window <= MAX_WINDOW) {
				// We've found a match, return
				//return true;
				char win[256];
				memset(win, 0, 256);
				strncpy(win, contig+*v+SEQ_LEN, window);
				printf("%s\n", win);
			}
		}
	}
}

void find_candidates(char* v_file, char* j_file, char* contig_file, int max_dist) {
	vmers = new sparse_hash_set<unsigned long, kmer_hash, eqkmer>();
	jmers = new sparse_hash_set<unsigned long, kmer_hash, eqkmer>();

	fprintf(stderr, "Loading vmers\n");
	fflush(stderr);
	load_kmers(v_file, vmers, max_dist);
	fprintf(stderr, "Loading jmers\n");
	fflush(stderr);
	load_kmers(j_file, jmers, max_dist);

//	char match = matches_vmer(5205);
//	printf("match1: %d\n", match);
//	match = matches_vmer(3000000000);
//	printf("match2: %d\n", match);

	fprintf(stderr, "Processing contigs\n");
	fflush(stderr);

	FILE* fp = fopen(contig_file, "r");

	char contig[100000];

	while (fgets(contig, 100000, fp) != NULL) {
		if(contig[0] != '>') {
			print_windows(contig);
		}
	}

	fclose(fp);

	fprintf(stderr, "Done");
	fflush(stderr);
}

int main(int argc, char** argv) {
	char* v_file = argv[1];
	char* j_file = argv[2];
	char* contig_file = argv[3];
	int max_dist = atoi(argv[4]);

	find_candidates(v_file, j_file, contig_file, max_dist);
}
