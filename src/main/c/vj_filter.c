#include "seq_dist.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <sparsehash/sparse_hash_map>
#include <sparsehash/sparse_hash_set>
#include <sparsehash/dense_hash_map>
#include <sparsehash/dense_hash_set>
#include "hash_utils.h"
#include "vj_filter.h"

using namespace std;
using google::sparse_hash_set;
using google::dense_hash_map;
using google::dense_hash_set;


int MIN_WINDOW;
int MAX_WINDOW;
int WINDOW_SPAN;
char J_CONSERVED;
int J_EXTENSION;

#define FRAME_PADDING 100
#define VJ_SEARCH_END 1000

// Anchors are located this many bases away from germline 3' and 5' ends of the V and J regions respectively
#define ANCHOR_PADDING 16

pthread_mutex_t vjf_mutex;

uint64_t MurmurHash64A ( const void * key, int len, uint64_t seed );

struct eqkmer
{
  bool operator()(unsigned long l1, unsigned long l2) const
  {
    return l1 == l2;
  }
};

struct kmer_hash
{
	uint64_t operator()(unsigned long kmer) const
	{
		return MurmurHash64A(&kmer, 4, 97);
		//return chunk;
	}
};

dense_hash_set<unsigned long, kmer_hash, eqkmer>* vmers;
dense_hash_set<unsigned long, kmer_hash, eqkmer>* jmers;

void load_kmers(char* input, dense_hash_set<unsigned long, kmer_hash, eqkmer>* kmers, int max_dist) {
	FILE* in = fopen(input, "r");

	unsigned long kmer;
	int freq;
	while (fscanf(in, "%lu\t%d\n", &kmer, &freq) == 2) {
		if (freq <= max_dist) {
			kmers->insert(kmer);
		}
	}

	fclose(in);
}

char matches_vmer(unsigned long kmer) {
	dense_hash_set<unsigned long, kmer_hash, eqkmer>::const_iterator it = vmers->find(kmer);
	return it != vmers->end();
}

char matches_jmer(unsigned long kmer) {
	dense_hash_set<unsigned long, kmer_hash, eqkmer>::const_iterator it = jmers->find(kmer);
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
		sparse_hash_set<const char*, vjf_hash, vjf_eqstr>& cdr3_seq, char*& cdr3) {
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
			fprintf(stderr, "C @ %d\n", i);
			v_indices.push_back(i);
		}
	}

	// TODO: Parameterize for Phenylalaline
	for (int i=j_index-ANCHOR_PADDING; i<SEQ_LEN+j_index-2; i++) {
		char ch[4];
		strncpy(ch, contig+i, 3);
		ch[3] = 0;

		if (J_CONSERVED == 'W' && strncmp(contig+i, "TGG", 3) == 0) {
			fprintf(stderr, "W @ %d\n", i);
			j_indices.push_back(i);
		} else if (J_CONSERVED == 'F' && (strncmp(contig+i, "TTC", 3) == 0 || strncmp(contig+i, "TTT", 3) == 0)) {
			fprintf(stderr, "J @ %d\n", i);
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

			fprintf(stderr, "window: %d\tj:\t%dv:%d\n", window, *j, *v);

			// TODO: Parameterize min/max window
			if (window % 3 == 0 && window >= MIN_WINDOW && window <= MAX_WINDOW && *j >= *v) {

				strncpy(cdr3, contig+*v, window);
				cdr3[window] = '\0';

				fprintf(stderr, "cdr3: %s\n", cdr3);

//				if (PRINT_CDR3_INDEX) {
//					fprintf(stderr, "contig: [%s]\t cdr3: [%s]\tV: [%d]\tJ: [%d]\n", contig, cdr3, *v, *j);
//					fprintf(stderr, "contig: [%s]\tV: [%d]\tJ: [%d]\n", contig, *v, *j);
//					printf("%s\t%s\t%d\t%d\n", contig, cdr3, *v, *j);
//				}

				cdr3_seq.insert(cdr3);
				cdr3 += 1024;
			}
		}
	}
}

char is_sub_string(const char* str, sparse_hash_set<const char*, vjf_hash, vjf_eqstr> cdr3_seq) {
	for (sparse_hash_set<const char*, vjf_hash, vjf_eqstr>::iterator it=cdr3_seq.begin(); it!=cdr3_seq.end(); it++) {
		if (strstr(*it, str) != NULL && strcmp(*it, str) !=0) {
			return true;
		}
	}

	return false;
}

// Thread local CDR3 buffer
__thread char* vjf_cdr3_block_buffer;

void print_windows(char* contig, dense_hash_map<const char*, const char*, vjf_hash, vjf_eqstr>& windows) {

	char* cdr3_block = vjf_cdr3_block_buffer;

	sparse_hash_set<const char*, vjf_hash, vjf_eqstr> cdr3_seq;

	char* contig_index = contig;
	vector<int> v_indices;
	vector<int> j_indices;

	int len = strlen(contig) - SEQ_LEN;

	for (int i=0; i<len; i++) {
		unsigned long kmer = seq_to_int(contig_index);

//		if (kmer == 0) {
//			printf("STRANGE CONTIG: %s\n", contig);
//			fflush(stdout);
//		}

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
			int window = *j - *v + SEQ_LEN;

			int pad = SEQ_LEN*2 + ANCHOR_PADDING*2;

			if (window >= MIN_WINDOW && window <= MAX_WINDOW+pad && strlen(contig) > window) {
				find_conserved_aminos(*v, *j, contig, cdr3_seq, cdr3_block);
			}
		}
	}

	sparse_hash_set<const char*, vjf_hash, vjf_eqstr>::const_iterator it;
	for (it=cdr3_seq.begin(); it!=cdr3_seq.end(); it++) {

		fprintf(stderr, "iterating: %s\n", *it);

		if (!is_sub_string(*it, cdr3_seq)) {

			int window = strlen(*it);
			fprintf(stderr, "window: %d\n", window);

			// Find current CDR3 string in contig
			char* start = strstr(contig, *it);

			if (start != NULL) {

				int vpad = WINDOW_SPAN - (window + J_EXTENSION);
				fprintf(stderr, "vpad: %d\n", vpad);

				start -= vpad;

				fprintf(stderr, "start length: %d\n", strlen(start));

				if (strlen(start) > WINDOW_SPAN) {
					char win[1024];
					memset(win, 0, 1024);
					strncpy(win, start, WINDOW_SPAN);

					fprintf(stderr, "Found win: %s\n", win);

					if (is_in_frame(win)) {
						fprintf(stderr, "Win in frame\n");
						if (windows.find(win) == windows.end()) {

							char* final_win = (char*) calloc(1024, sizeof(char));
							char* final_cdr3 = (char*) calloc(256, sizeof(char));
							strcpy(final_win, win);
							strcpy(final_cdr3, *it);

							fprintf(stderr, "Adding to windows: %s\t%s\n", final_win, final_cdr3);

							windows[final_win] = final_cdr3;

							fprintf(stderr, "windows size: %d\n", windows.size());
						}
					}
				}
			}
		}
	}
}

void init(char* v_file, char* j_file, int max_dist) {
	pthread_mutex_init(&vjf_mutex, NULL);

	vmers = new dense_hash_set<unsigned long, kmer_hash, eqkmer>();
	jmers = new dense_hash_set<unsigned long, kmer_hash, eqkmer>();

	vmers->set_empty_key(0);
	jmers->set_empty_key(0);

	fprintf(stderr, "Loading vmers\n");
	fflush(stderr);
	load_kmers(v_file, vmers, max_dist);
	fprintf(stderr, "Loading jmers\n");
	fflush(stderr);
	load_kmers(j_file, jmers, max_dist);

}

void vjf_init(char* v_file, char* j_file, int max_dist, int min_win, int max_win,
		char j_conserved, int window_span, int j_extension) {

	MIN_WINDOW = min_win;
	MAX_WINDOW = max_win;
	J_CONSERVED = j_conserved;
	WINDOW_SPAN = window_span;
	J_EXTENSION = j_extension;

	init(v_file, j_file, max_dist);
}



void vjf_search(char* contig, dense_hash_map<const char*, const char*, vjf_hash, vjf_eqstr>& windows) {
	print_windows(contig, windows);
}

/*
void find_candidates(char* v_file, char* j_file, char* contig_file, int max_dist) {

	init(v_file, j_file, max_dist);

	fprintf(stderr, "Processing contigs\n");
	fflush(stderr);

	FILE* fp = fopen(contig_file, "r");

	char contig[100000];

	dense_hash_set<const char*, vjf_hash, vjf_eqstr> windows;

	int i = 0;
	while (fgets(contig, 100000, fp) != NULL) {
		if(contig[0] != '>') {
//			memset(cdr3_block, 0, 256*1000000);
			vjf_search(contig, windows);
			if (i++ % 10000 == 0) {
				fprintf(stderr, "Processed %d contigs\n", i);
			}
		}
	}

	fprintf(stderr, "Outputting windows\n");

	for (dense_hash_set<const char*, vjf_hash, vjf_eqstr>::iterator it=windows.begin(); it!=windows.end(); it++) {
		printf("%s\n", *it);
	}

	fclose(fp);

	fprintf(stderr, "Done");
	fflush(stderr);
}
*/

/*
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
	J_EXTENSION = atoi(argv[9]);

	fprintf(stderr, "V: %s\nJ: %s\ncontigs: %s\ndist: %d\nmin_win: %d\nmax_win: %d\nj_cons: %c\nwindow_span: %d\nj_extension: %d\n",
			v_file, j_file, contig_file, max_dist, MIN_WINDOW, MAX_WINDOW, J_CONSERVED, WINDOW_SPAN, J_EXTENSION);

	find_candidates(v_file, j_file, contig_file, max_dist);
}
*/


