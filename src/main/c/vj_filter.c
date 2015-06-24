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
//#define MIN_WINDOW 10
//#define MAX_WINDOW 90

//#define WINDOW_SPAN 90

int MIN_WINDOW;
int MAX_WINDOW;
int WINDOW_SPAN;
char J_CONSERVED;

#define FRAME_PADDING 100
#define VJ_SEARCH_END 1000

// Anchors are located this many bases away from germline 3' and 5' ends of the V and J regions respectively
#define ANCHOR_PADDING 10

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

struct eqstr
{
  bool operator()(const char* s1, const char* s2) const
  {
    return (s1 == s2) || (s1 && s2 && strcmp(s1, s2) == 0);
  }
};

struct my_hash
{
	uint64_t operator()(const char* seq) const
	{
		uint64_t h = MurmurHash64A(seq, strlen(seq), 97);
		return h;
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

char is_stop_codon(char* codon) {
	return strncmp(codon, "TAG", 3) == 0 || strncmp(codon, "TAA", 3) == 0 || strncmp(codon, "TGA", 3) == 0;
}

char is_in_frame(int start, int end, int frame, char* contig) {
	for (int i=start+frame; i<end+frame-3; i+=3) {
		char* codon = contig + i;
		if (is_stop_codon(codon)) {
			return false;
		}
	}

	return true;
}

char is_in_frame(char* seq) {
	for (int i=0; i<strlen(seq)-2; i+=3) {
		char* codon = seq + i;
		if (is_stop_codon(codon)) {
			return false;
		}
	}

	return true;
}

// Look for a window of a few hundred bases without a stop codon
int find_frame(char* contig, int v, int window) {

	int start = v - FRAME_PADDING;
	int end = v + window + FRAME_PADDING;

	int frame = -1;;

	for (int i=0; i<3; i++) {
		if (is_in_frame(start, end, i, contig)) {
			if (frame != -1) {
				fprintf(stderr, "Multiple frames in: %s : %d\n", contig, v);
			} else {
				frame = i;
			}
		}
	}

	return frame;
}

void find_conserved_aminos(int v_index, int j_index, char* contig,
		sparse_hash_set<const char*, my_hash, eqstr>& cdr3_seq, char*& cdr3) {
	// BCR heavy has conserved Cysteine in V and Tryptophan in J
	// BCR light, TCR light & heavy have conserved Cysteine in V and Phenylalaline in J
	// Cysteine : C :  TGT,TGC
	// Tryptophan : W : TGG
	// Phenylalaline : F : TTC,TTT
	//
	// See: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC441550/

	vector<int> v_indices;
	vector<int> j_indices;

	for (int i=v_index; i<SEQ_LEN+ANCHOR_PADDING+v_index-2; i++) {
		char ch[4];
		strncpy(ch, contig+i, 3);
		ch[3] = 0;

		if (strncmp(contig+i, "TGT", 3) == 0 || strncmp(contig+i, "TGC", 3) == 0) {
//			fprintf(stderr, "C @ %d\n", i);
			v_indices.push_back(i);
		}
	}

	// TODO: Parameterize for Phenylalaline
	for (int i=j_index-ANCHOR_PADDING; i<SEQ_LEN+j_index-2; i++) {
		char ch[4];
		strncpy(ch, contig+i, 3);
		ch[3] = 0;

		if (J_CONSERVED == 'W' && strncmp(contig+i, "TGG", 3) == 0) {
//			fprintf(stderr, "W @ %d\n", i);
			j_indices.push_back(i);
		} else if (J_CONSERVED == 'F' && (strncmp(contig+i, "TTC", 3) == 0 || strncmp(contig+i, "TTT", 3) == 0)) {
//			fprintf(stderr, "J @ %d\n", i);
			j_indices.push_back(i);
		}
	}

//	char* cdr3 = cdr3_block;

	vector<int>::const_iterator v;
	for (v=v_indices.begin(); v!=v_indices.end(); v++) {
		vector<int>::const_iterator j;
		for (j=j_indices.begin(); j!=j_indices.end(); j++) {
			// For each C, (W/F) combo, check the distance and frame
			int window = *j - *v + 3;

			// TODO: Parameterize min/max window
			if (window % 3 == 0 && window >= MIN_WINDOW && window <= MAX_WINDOW && *j >= *v) {

				strncpy(cdr3, contig+*v, window);
				cdr3_seq.insert(cdr3);
				cdr3 += 256;
/*

				int pad = MAX_WINDOW - window;
				int vpad = pad / 2;
				vpad -= vpad % 3;
				int start = *v - vpad;

				if (strlen(contig) > start+MAX_WINDOW) {
					char win[256];
					memset(win, 0, 256);
					strncpy(win, contig+start, MAX_WINDOW);
					printf("%s\n", win);
				}
*/
			}

		}
	}
}

char is_sub_string(const char* str, sparse_hash_set<const char*, my_hash, eqstr> cdr3_seq) {
	for (sparse_hash_set<const char*, my_hash, eqstr>::iterator it=cdr3_seq.begin(); it!=cdr3_seq.end(); it++) {
		if (strstr(*it, str) != NULL && strcmp(*it, str) !=0) {
			return true;
		}
	}

	return false;
}

void print_windows(char* contig, sparse_hash_set<const char*, my_hash, eqstr>& windows) {

	// Allocate space for up to 1,000,000 CDR3's
	char* cdr3_block = (char*) calloc(256*1000000, sizeof(char));
	char* orig_cdr3_block = cdr3_block;

	sparse_hash_set<const char*, my_hash, eqstr> cdr3_seq;

	char* contig_index = contig;
	vector<int> v_indices;
	vector<int> j_indices;

	int len = strlen(contig) - SEQ_LEN;

	for (int i=0; i<len; i++) {
		unsigned long kmer = seq_to_int(contig_index);

		if (matches_vmer(kmer)) {
//			fprintf(stderr, "v: %d\n", i);
			v_indices.push_back(i);
		}

		if (matches_jmer(kmer)) {
//			fprintf(stderr, "j: %d\n", i);
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
			int window = *j - *v + SEQ_LEN;

			int pad = SEQ_LEN*2 + ANCHOR_PADDING*2;

			if (window >= MIN_WINDOW && window <= MAX_WINDOW+pad && strlen(contig) > window) {
				find_conserved_aminos(*v, *j, contig, cdr3_seq, cdr3_block);
				// We've found a potential match
//				char win[256];
//				memset(win, 0, 256);
//				strncpy(win, contig+*v, window);
//				printf("%s\n", win);
			}
		}
	}

	sparse_hash_set<const char*, my_hash, eqstr>::const_iterator it;
	for (it=cdr3_seq.begin(); it!=cdr3_seq.end(); it++) {

		if (!is_sub_string(*it, cdr3_seq)) {

			int window = strlen(*it);

			// Find current CDR3 string in contig
			char* start = strstr(contig, *it);

			if (start != NULL) {
				int pad = WINDOW_SPAN - window;
				int vpad = pad / 2;
				if (vpad % 3 == 1) {
					vpad += 2;
				} else if (vpad % 3 == 2) {
					vpad += 1;
				}
				start -= vpad;

				/*
				int pad = MAX_WINDOW - window;
				int vpad = pad / 2;
				vpad -= vpad % 3;
				start -= vpad;
				*/

				if (strlen(start) > WINDOW_SPAN) {
					char win[256];
					memset(win, 0, 256);
					strncpy(win, start, WINDOW_SPAN);

					if (is_in_frame(win)) {
						if (windows.find(win) == windows.end()) {

							char* final_win = (char*) calloc(256, sizeof(char));
							strcpy(final_win, win);

							windows.insert(final_win);
						}
					}
				}
			}
		}
	}

	free(orig_cdr3_block);
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

	sparse_hash_set<const char*, my_hash, eqstr> windows;

	int i = 0;
	while (fgets(contig, 100000, fp) != NULL) {
		if(contig[0] != '>') {
//			memset(cdr3_block, 0, 256*1000000);
			print_windows(contig, windows);
			if (i++ % 10000 == 0) {
				fprintf(stderr, "Processed %d contigs\n", i);
			}
		}
	}

	fprintf(stderr, "Outputting windows\n");

	for (sparse_hash_set<const char*, my_hash, eqstr>::iterator it=windows.begin(); it!=windows.end(); it++) {
		printf("%s\n", *it);
	}

	fclose(fp);

	fprintf(stderr, "Done");
	fflush(stderr);
}

int main(int argc, char** argv) {
	fprintf(stderr, "Here...\n");

	char* v_file = argv[1];
	char* j_file = argv[2];
	char* contig_file = argv[3];
	int max_dist = atoi(argv[4]);
	MIN_WINDOW = atoi(argv[5]);
	MAX_WINDOW = atoi(argv[6]);
	J_CONSERVED = argv[7][0];
	WINDOW_SPAN = atoi(argv[8]);

	fprintf(stderr, "V: %s\nJ: %s\ncontigs: %s\ndist: %d\nmin_win: %d\nmax_win: %d\nj_cons: %c\nwindow_span: %d\n",
			v_file, j_file, contig_file, max_dist, MIN_WINDOW, MAX_WINDOW, J_CONSERVED, WINDOW_SPAN);

	find_candidates(v_file, j_file, contig_file, max_dist);
}

