/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <pthread.h>
#include <unistd.h>
#include <iostream>
#include <stack>
#include <list>
#include <vector>
#include <sparsehash/sparse_hash_map>
#include <sparsehash/sparse_hash_set>
#include <stdexcept>
//#include "abra_NativeAssembler.h"
#include "seq_score.h"

using namespace std;
using google::sparse_hash_map;
using google::sparse_hash_set;

//#define READ_LENGTH 100
//#define KMER 63
//#define MIN_CONTIG_LENGTH 101
//#define MIN_NODE_FREQUENCY 3
//#define MIN_NODE_FREQUENCY 2
#define MAX_CONTIG_SIZE 10000
#define MAX_READ_LENGTH 1001
//#define MIN_BASE_QUALITY 20
#define INCREASE_MIN_NODE_FREQ_THRESHOLD 2500

#define MAX_TOTAL_CONTIG_LEN 10000000

#define OK 0
#define TOO_MANY_PATHS_FROM_ROOT -1
#define TOO_MANY_CONTIGS -2
#define STOPPED_ON_REPEAT -3
#define TOO_MANY_NODES -4

#define MAX_FREQUENCY 32766
#define MAX_QUAL_SUM 255

//#define MIN_HOMOLOGY_SCORE 13

// TODO: This is used to bound qual sum arrays.  Use a memory pool instead for this.
#define MAX_KMER_LEN 201

// This makes sense for small assembly windows, but should be parameterized for larger assemblies
#define MAX_NODES 900000000

// Kmers containing bases below this threshold are excluded from assembly.
#define MIN_BASE_QUALITY 13

// Minimum edge frequency as percent
// #define MIN_EDGE_FREQUENCY .02

#define MIN_EDGE_FREQUENCY -1

// TODO: This should vary be kmer len
#define MIN_ROOT_HOMOLOGY_SCORE 16

// Must be greater than the number of  source(root) nodes - TODO: re-use threads and allocate dynamically.
#define MAX_THREADS 100000

// TODO: Parameterize
int MAX_RUNNING_THREADS;

pthread_mutex_t running_thread_mutex;
pthread_mutex_t contig_writer_mutex;
pthread_mutex_t marker_trackback_mutex;

int running_threads = 0;

#define MAX_DISTANCE_FROM_MARKER 1000

// #define MAX_TRACEBACK_STACK_SIZE 1024
#define MAX_TRACEBACK_STACK_SIZE 256


//TODO: Better variable localization
//__thread int read_length;
//__thread int min_contig_length;
//__thread int kmer_size;
//__thread int min_node_freq;
//__thread int min_base_quality;

int read_length;
int min_contig_length;
int kmer_size;
int min_node_freq;
int min_base_quality;

//#define MIN_CONTIG_SCORE -6
float MIN_CONTIG_SCORE;


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


struct eqstr
{
  bool operator()(const char* s1, const char* s2) const
  {
    return (s1 == s2) || (s1 && s2 && strncmp(s1, s2, kmer_size) == 0);
  }
};

struct my_hash
{
	uint64_t operator()(const char* kmer) const
	{
		return MurmurHash64A(kmer, kmer_size, 97);
		//return chunk;
	}
};

struct struct_pool {
	struct node_pool* node_pool;
	struct read_pool* read_pool;
};

#define NODES_PER_BLOCK 1000000
#define MAX_NODE_BLOCKS 50000
#define READS_PER_BLOCK 1000000
#define MAX_READ_BLOCKS 100000

struct node_pool {
	struct node** nodes;
	int block_idx;
	int node_idx;
};

struct read_pool {
	char** reads;
	int block_idx;
	int read_idx;
};

struct node {

	//TODO: Collapse from 8 to 2 bits.  Only store as key.
	char* kmer;
	//TODO: Convert to stl?
	struct linked_node* toNodes;
	struct linked_node* fromNodes;
	char* contributingRead;
	unsigned char qual_sums[MAX_KMER_LEN];
	unsigned short frequency;
	char hasMultipleUniqueReads;
	char contributing_strand;
	char root_eligible;
	char is_simplified;
	char is_root;
};

struct linked_node {
	struct node* node;
	struct linked_node* next;
};

int compare(const char* s1, const char* s2) {
	return (s1 == s2) || (s1 && s2 && strcmp(s1, s2) == 0);
}

int compare_kmer(const char* s1, const char* s2) {
	return (s1 == s2) || (s1 && s2 && strncmp(s1, s2, kmer_size) == 0);
}

struct struct_pool* init_pool() {

	printf("Memory pool init start\n");
	fflush(stdout);
	struct_pool* pool = (struct struct_pool*) malloc(sizeof(struct_pool));
	pool->node_pool = (struct node_pool*) malloc(sizeof(node_pool));
	// Allocate array of arrays
	pool->node_pool->nodes = (struct node**) malloc(sizeof(struct node*) * MAX_NODE_BLOCKS);
	// Allocate first array of nodes
	pool->node_pool->nodes[0] = (struct node*) malloc(sizeof(struct node) * NODES_PER_BLOCK);
	pool->node_pool->block_idx = 0;
	pool->node_pool->node_idx = 0;

	pool->read_pool = (struct read_pool*) malloc(sizeof(read_pool));
	pool->read_pool->reads = (char**) malloc(sizeof(char*) * MAX_READ_BLOCKS);
	pool->read_pool->reads[0] = (char*) malloc(sizeof(char) * (read_length+1) * READS_PER_BLOCK);
	pool->read_pool->block_idx = 0;
	pool->read_pool->read_idx = 0;
	printf("Memory pool init done\n");
	fflush(stdout);

	return pool;
}

char* allocate_read(struct_pool* pool) {
	if (pool->read_pool->block_idx > MAX_READ_BLOCKS) {
		printf("READ BLOCK INDEX TOO BIG!!!!\n");
		exit(-1);
	}

	if (pool->read_pool->read_idx >= READS_PER_BLOCK) {
		pool->read_pool->block_idx++;
		pool->read_pool->read_idx = 0;
		pool->read_pool->reads[pool->read_pool->block_idx] = (char*) malloc(sizeof(char) * (read_length+1) * READS_PER_BLOCK);
	}

	return &pool->read_pool->reads[pool->read_pool->block_idx][pool->read_pool->read_idx++ * (read_length+1)];
}

struct node* allocate_node(struct_pool* pool) {
	if (pool->node_pool->block_idx >= MAX_NODE_BLOCKS) {
		printf("NODE BLOCK INDEX TOO BIG!!!!\n");
		exit(-1);
	}

	if (pool->node_pool->node_idx >= NODES_PER_BLOCK) {
printf("Increasing node pool...\n");
fflush(stdout);
		pool->node_pool->block_idx++;
		pool->node_pool->node_idx = 0;
		pool->node_pool->nodes[pool->node_pool->block_idx] = (struct node*) malloc(sizeof(struct node) * NODES_PER_BLOCK);
	}

	return &pool->node_pool->nodes[pool->node_pool->block_idx][pool->node_pool->node_idx++];
}

unsigned char phred33(char ch) {
	return ch - '!';
}

struct node* new_node(char* seq, char* contributingRead, struct_pool* pool, int strand, char* quals) {

//	node* my_node = (node*) malloc(sizeof(node));
	node* my_node = allocate_node(pool);
	memset(my_node, 0, sizeof(node));
	my_node->kmer = seq;
//	strcpy(my_node->contributingRead, contributingRead);
	my_node->contributingRead = contributingRead;
	my_node->frequency = 1;
	my_node->hasMultipleUniqueReads = 0;
	my_node->contributing_strand = (char) strand;
	for (int i=0; i<kmer_size; i++) {
		my_node->qual_sums[i] = phred33(quals[i]);
	}
	return my_node;
}

char* get_kmer(int idx, char* sequence) {
	return &sequence[idx];
}

int is_node_in_list(struct node* node, struct linked_node* list) {
	struct linked_node* ptr = list;

	while (ptr != NULL) {
		if (compare_kmer(ptr->node->kmer, node->kmer)) {
			return 1;
		}
		ptr = ptr->next;
	}

	return 0;
}

void link_nodes(struct node* from_node, struct node* to_node) {
	if (!is_node_in_list(to_node, from_node->toNodes)) {
		struct linked_node* to_link = (linked_node*) malloc(sizeof(linked_node));
		to_link->node = to_node;
		to_link->next = from_node->toNodes;
		from_node->toNodes = to_link;
	}

	if (!is_node_in_list(from_node, to_node->fromNodes)) {
		struct linked_node* from_link = (linked_node*) malloc(sizeof(linked_node));
		from_link->node = from_node;
		from_link->next = to_node->fromNodes;
		to_node->fromNodes = from_link;
	}
}

void increment_node_freq(struct node* node, char* read_seq, int strand, char* kmer_qual) {
	if (node->frequency < MAX_FREQUENCY-1) {
		node->frequency++;
	} else {
//		printf("Max frequency reached for: %s\n", node->kmer);
	}

	if (!(node->hasMultipleUniqueReads) &&
		(!compare(node->contributingRead, read_seq) || node->contributing_strand != (char) strand)) {
		node->hasMultipleUniqueReads = 1;
	}

	for (int i=0; i<kmer_size; i++) {
		unsigned char phred33_qual = phred33(kmer_qual[i]);
		if ((node->qual_sums[i] + phred33_qual) < MAX_QUAL_SUM) {
			node->qual_sums[i] += phred33_qual;
		} else {
			node->qual_sums[i] = MAX_QUAL_SUM;
		}
	}
}

int include_kmer(char* sequence, char*qual, int idx) {

	int include = 1;

	for (int i=idx; i<idx+kmer_size; i++) {
		// Discard kmers with ambiguous bases
		if (sequence[i] == 'N') {
			include = 0;
			break;
		}

		// Discard kmers with low base qualities

//		if (qual[i] - '!' < min_base_quality) {
		if (phred33(qual[i]) < MIN_BASE_QUALITY) {
			include = 0;
			break;
		}
	}

	return include;
}

void add_to_graph(char* sequence, sparse_hash_map<const char*, struct node*, my_hash, eqstr>* nodes, struct_pool* pool, char* qual, int strand, char has_roots) {

	struct node* prev = 0;

	for (int i=0; i<=read_length-kmer_size; i++) {

		if (include_kmer(sequence, qual, i)) {
			char* kmer = get_kmer(i, sequence);
			char* kmer_qual = get_kmer(i, qual);

			struct node* curr = (*nodes)[kmer];

			if (curr == NULL) {
				curr = new_node(kmer, sequence, pool, strand, kmer_qual);
				curr->root_eligible = has_roots;

				if (curr == NULL) {
					printf("Null node for kmer: %s\n", kmer);
					exit(-1);
				}

				(*nodes)[kmer] = curr;
			} else {
				if (has_roots && !curr->root_eligible) {
					curr->root_eligible = 1;
				}
				increment_node_freq(curr, sequence, strand, kmer_qual);
			}

			if (prev != NULL) {
				link_nodes(prev, curr);
			}

			prev = curr;
		} else {
			prev = NULL;
		}
	}
}

void build_graph2(const char* input, sparse_hash_map<const char*, struct node*, my_hash, eqstr>* nodes, struct_pool* pool, char has_roots) {
	size_t input_len = strlen(input);
	int record_len = read_length*2 + 1;

printf("input_len: %ld, record_len: %d\n", input_len, record_len);
	size_t num_records = input_len / record_len;
	size_t record = 0;
	const char* ptr = input;
	int num_reads = 0;

	while ((record < num_records) && (nodes->size() < MAX_NODES)) {
		ptr = &(input[record*record_len]);
		int strand = 0;

		if (ptr[0] == '0') {
			strand = 0;
		} else if (ptr[0] == '1') {
			strand = 1;
		} else {
			printf("Initial char in input invalid: %c\n", ptr[0]);
			printf("ERROR!  INVALID INPUT:\n===========================%s\n===========================\n", input);
			exit(-1);
		}

		// TODO: skip copying the input read.  Downstream code appears to depend
		// upon null terminator.
		char* read_ptr = allocate_read(pool);
		memset(read_ptr, 0, read_length+1);
		memcpy(read_ptr, &(ptr[1]), read_length);

//		char* read_ptr = (char*) &(ptr[1]);

		char* qual_ptr = (char*) &(ptr[read_length+1]);
		add_to_graph(read_ptr, nodes, pool, qual_ptr, strand, has_roots);
		record++;

		if ((record % 1000000) == 0) {
			printf("record_count: %d\n", record);
			fflush(stdout);
		}
	}

	printf("Num reads: %d\n", record);
	printf("Num nodes: %d\n", nodes->size());
	fflush(stdout);
}

/*
void build_graph(const char* read_file, sparse_hash_map<const char*, struct node*, my_hash, eqstr>* nodes, struct_pool* pool) {
	FILE *fp = fopen(read_file, "r");
	char read[MAX_READ_LENGTH];
	memset(read, 0, MAX_READ_LENGTH);

	char qual[MAX_READ_LENGTH];
	memset(qual, 0, MAX_READ_LENGTH);

	char read_info[128];

	int line = 0;
	while (fscanf(fp, "%s", read_info) != EOF) {
		fscanf(fp, "%s", read);
		fscanf(fp, "%s", qual);
		if (strcmp(read, "") != 0) {
			char* read_ptr = allocate_read(pool);
			memcpy(read_ptr, read, read_length+1);
			int strand = atoi(read_info);
			add_to_graph(read_ptr, nodes, pool, qual, strand);
			line++;

			if ((line % 100000) == 0) {
				printf("Processed %d reads.\n", line);
			}
		}
	}

	printf("Num reads: %d\n", line);
	printf("Num nodes: %d\n", nodes->size());

	fclose(fp);
}
*/

struct linked_node* remove_node_from_list(struct node* node, struct linked_node* list) {
	struct linked_node* node_ptr = list;
	struct linked_node* prev_ptr = NULL;

	char is_found = false;
	while ((node_ptr != NULL) && (!is_found)) {
		if (strncmp(node_ptr->node->kmer, node->kmer, kmer_size) == 0) {
			if (prev_ptr == NULL) {
				// Set head of list to next elem
				list = list->next;
			} else {
				// Remove node from list
				prev_ptr->next = node_ptr->next;
			}

			// Free linked_node
			free(node_ptr);
			is_found = true;
		}

		prev_ptr = node_ptr;
		node_ptr = node_ptr->next;
	}

	return list;
}

void cleanup(struct linked_node* linked_nodes) {
	struct linked_node* ptr = linked_nodes;
	while (ptr != NULL) {
		struct linked_node* next = ptr->next;
		free(ptr);
		ptr = next;
	}
}

int is_base_quality_good(struct node* node) {
	int is_good = 1;

	for (int i=0; i<kmer_size; i++) {
		if (node->qual_sums[i] < min_base_quality) {
			is_good = 0;
			break;
		}
	}

	return is_good;
}

void remove_node_and_cleanup(const char* key, struct node* node, sparse_hash_map<const char*, struct node*, my_hash, eqstr>* nodes) {

	if (!node->is_root) {
		// Remove node from "from" lists
		struct linked_node* to_node = node->toNodes;
		while (to_node != NULL) {
			to_node->node->fromNodes = remove_node_from_list(node, to_node->node->fromNodes);
			to_node = to_node->next;
		}

		// Remove node from "to" lists
		struct linked_node* from_node = node->fromNodes;
		while (from_node != NULL) {
			from_node->node->toNodes = remove_node_from_list(node, from_node->node->toNodes);
			from_node = from_node->next;
		}

		// Remove node from map
		nodes->erase(key);
		cleanup(node->toNodes);
		node->toNodes = NULL;
		cleanup(node->fromNodes);
		node->fromNodes = NULL;
	}
}

void prune_low_frequency_edges(sparse_hash_map<const char*, struct node*, my_hash, eqstr>* nodes) {

	long removed_edge_count = 0;

	for (sparse_hash_map<const char*, struct node*, my_hash, eqstr>::const_iterator it = nodes->begin();
				 it != nodes->end(); ++it) {

		const char* key = it->first;
		node* curr_node = it->second;

		if (curr_node != NULL) {
			////////////////////////////////////////////////
			// Check to node list for low frequency edges
			struct linked_node* to_node = curr_node->toNodes;

			// Calculate total outgoing "edge" frequency
			int to_node_total_freq = 0;

			while (to_node != NULL) {
				// Using node frequency as proxy for edge frequency here...
				to_node_total_freq = to_node_total_freq + to_node->node->frequency;
				to_node = to_node->next;
			}

			// Identify edges to prune
			to_node = curr_node->toNodes;
			vector<node*> to_nodes_to_remove;

			while (to_node != NULL) {
				if ( ((double) to_node->node->frequency / (double) to_node_total_freq) < MIN_EDGE_FREQUENCY ) {
					to_nodes_to_remove.push_back(to_node->node);
				}
				to_node = to_node->next;
			}

			// Remove edges
			for (vector<node*>::const_iterator iter = to_nodes_to_remove.begin(); iter != to_nodes_to_remove.end(); ++iter ) {
				// Remove edges in each direction
				node* node_to_remove = *(iter);
				node_to_remove->fromNodes = remove_node_from_list(curr_node, node_to_remove->fromNodes);
				curr_node->toNodes = remove_node_from_list(node_to_remove, curr_node->toNodes);
				removed_edge_count += 1;
			}

			////////////////////////////////////////////////
			// Check from node list for low frequency edges
			struct linked_node* from_node = curr_node->fromNodes;

			// Calculate total outgoing "edge" frequency
			int from_node_total_freq = 0;

			while (from_node != NULL) {
				// Using node frequency as proxy for edge frequency here...
				from_node_total_freq = from_node_total_freq + from_node->node->frequency;
				from_node = from_node->next;
			}

			// Identify edges to prune
			from_node = curr_node->fromNodes;
			vector<node*> from_nodes_to_remove;

			while (from_node != NULL) {
				if ( ((double) from_node->node->frequency / (double) from_node_total_freq) < MIN_EDGE_FREQUENCY ) {
					from_nodes_to_remove.push_back(from_node->node);
				}
				from_node = from_node->next;
			}

			// Remove edges
			for (vector<node*>::const_iterator iter = from_nodes_to_remove.begin(); iter != from_nodes_to_remove.end(); ++iter ) {
				// Remove edges in each direction
				node* node_to_remove = *(iter);
				node_to_remove->toNodes = remove_node_from_list(curr_node, node_to_remove->toNodes);
				curr_node->fromNodes = remove_node_from_list(node_to_remove, curr_node->fromNodes);
				removed_edge_count += 1;
			}
		}
	}

	printf("Pruned %ld edges\n", removed_edge_count);
}

void prune_graph(sparse_hash_map<const char*, struct node*, my_hash, eqstr>* nodes, char isUnalignedRegion) {

	// First prune kmers that do not reach base quality sum threshold
	for (sparse_hash_map<const char*, struct node*, my_hash, eqstr>::const_iterator it = nodes->begin();
				 it != nodes->end(); ++it) {

		const char* key = it->first;
		struct node* node = it->second;

		if (node != NULL && !is_base_quality_good(node)) {
			remove_node_and_cleanup(key, node, nodes);
		}
	}

	printf("Remaining nodes after pruning step 1: %d\n", nodes->size());

	// Now go back through and ensure that each node reaches minimum frequency threshold.
	int freq = min_node_freq;

	if (!isUnalignedRegion) {
		//int increase_freq = nodes->size() / INCREASE_MIN_NODE_FREQ_THRESHOLD;
		int increase_freq = 0;

		if (increase_freq > 0) {
			freq = freq + increase_freq;
			printf("Increased mnf to: %d for nodes size: %d\n", freq, nodes->size());
		}
	}

	if (freq > 1) {
                printf("Pruning nodes < frequency: %d\n", freq);
		for (sparse_hash_map<const char*, struct node*, my_hash, eqstr>::const_iterator it = nodes->begin();
					 it != nodes->end(); ++it) {

			const char* key = it->first;
			struct node* node = it->second;

			if ((node != NULL) && ((node->frequency < freq) || (!(node->hasMultipleUniqueReads)))) {
				remove_node_and_cleanup(key, node, nodes);
			}
		}
	}

	printf("Remaining nodes after pruning step 2: %d\n", nodes->size());

	prune_low_frequency_edges(nodes);

	// Final pass through cleaning up nodes that are unreachable
	for (sparse_hash_map<const char*, struct node*, my_hash, eqstr>::const_iterator it = nodes->begin();
				 it != nodes->end(); ++it) {

		const char* key = it->first;
		struct node* node = it->second;

		if (node != NULL && node->toNodes == NULL && node->fromNodes == NULL) {
			remove_node_and_cleanup(key, node, nodes);
		}
	}

	printf("Remaining nodes after edge pruning: %d\n", nodes->size());

}

void print_kmer(char* kmer) {
    for (int i=0; i<kmer_size; i++) {
            printf("%c", kmer[i]);
    }
}

void print_kmer(struct node* node) {
	printf("%x\t", node);
    for (int i=0; i<kmer_size; i++) {
            printf("%c", node->kmer[i]);
    }
}

void print_node(struct node* node) {
        printf("kmer: ");
        print_kmer(node);
        printf("\tfrom: ");

        struct linked_node* from = node->fromNodes;
        while (from != NULL) {
        	print_kmer(from->node);
        	printf(",");
        	from = from->next;
        }

        printf("\tto: ");
        struct linked_node* to = node->toNodes;
        while (to != NULL) {
        	print_kmer(to->node);
        	printf(",");
        	to = to->next;
        }
}

//#define SOURCE "AAACGGGCCGTTTGCATTGTGAACT"
#define SOURCE "GGGGCAGCAAGATGGTGTTGCAGACCCAGGT"

int num_root_candidates = 0;

int is_root(struct node* node, int& num_root_candidates) {
	int is_root = 0;

//	return node != NULL && strncmp(node->kmer, SOURCE, kmer_size) == 0;

	if (node != NULL && node->root_eligible) {
		if (node->fromNodes == NULL) {
			num_root_candidates += 1;
/*
			int homology_score = score_seq(node->kmer);
			printf("HOMOLOGY_SCORE:\t%d\t%d\t", homology_score, node->frequency);
			print_kmer(node->kmer);
			printf("\n");
			if (homology_score >= MIN_ROOT_HOMOLOGY_SCORE) {
				is_root = 1;
			}
*/
			if ((node->frequency >= min_node_freq) && (node->hasMultipleUniqueReads) && is_base_quality_good(node)) {
//			if ((node->frequency >= min_node_freq) && (node->hasMultipleUniqueReads) && is_base_quality_good(node) &&
//					(strncmp(node->kmer, "CGCACCATCTCCAAGGACACCTCCA", kmer_size) == 0) ) {

				is_root = 1;
				node->is_root = 1;

				printf("ROOT_NODE:\t%d\t", node->frequency);
				print_kmer(node->kmer);
				printf("\n");
			}
		} else {
			// Identify nodes that point to themselves with no other incoming edges.
			// This will be cleaned up during contig building.
			struct linked_node* from = node->fromNodes;
			if (from->next == NULL && (strncmp(node->kmer, from->node->kmer, kmer_size) == 0)) {
//				is_root = 1;
				printf("SELF_ROOT\n");
			}
		}
	}

	return is_root;

}

struct linked_node* identify_root_nodes(sparse_hash_map<const char*, struct node*, my_hash, eqstr>* nodes) {

	struct linked_node* root_nodes = NULL;
	int count = 0;
	int num_root_candidates;

	for (sparse_hash_map<const char*, struct node*, my_hash, eqstr>::const_iterator it = nodes->begin();
	         it != nodes->end(); ++it) {
		struct node* node = it->second;

//		if (node != NULL) {
//			print_node(node);
//		}

		if (is_root(node, num_root_candidates)) {
			struct linked_node* next = root_nodes;
			root_nodes = (linked_node*) malloc(sizeof(linked_node));
			root_nodes->node = node;
			root_nodes->next = next;

//			printf("\tROOT");

			count++;
		}

//		printf("\n");
	}

	printf("num root nodes: %d\n", count);

	printf("NUM_ROOT_CANDIDATES: %d\n", num_root_candidates);
	fflush(stdout);

	return root_nodes;
}

struct node_tracker {
	struct node* curr_node;
	sparse_hash_set<const char*, my_hash, eqstr>* visited_nodes;
	int distance;
};

void get_roots_from_marker(sparse_hash_map<const char*, struct node*, my_hash, eqstr>& roots, struct node* marker) {

	node_tracker* tracker = (node_tracker*) calloc(1, sizeof(node_tracker));
	tracker->curr_node = marker;
	tracker->visited_nodes = new sparse_hash_set<const char*, my_hash, eqstr>();
	tracker->distance = 0;
//	tracker->visited_nodes->insert(tracker.curr_node->kmer);

	stack<node_tracker*> node_stack;
	node_stack.push(tracker);

	printf("Track back start: ");
	print_kmer(tracker->curr_node);
	printf("\n");
	fflush(stdout);

	while (!node_stack.empty()) {
		node_tracker* curr = node_stack.top();
		node_stack.pop();

		// Search for the current node's kmer in the set of visited kmers
		sparse_hash_set<const char*, my_hash, eqstr>::const_iterator it = curr->visited_nodes->find(curr->curr_node->kmer);
		char is_repeat = it != curr->visited_nodes->end();

		if (node_stack.size() > MAX_TRACEBACK_STACK_SIZE) {
			printf("Traceback stack too big, skipping: ");
			print_kmer(tracker->curr_node);
			printf("\n");
			fflush(stdout);
			// Stop processing and clear stack
			delete curr->visited_nodes;
			free(curr);
			while (!node_stack.empty()) {
				node_tracker* curr = node_stack.top();
				node_stack.pop();
				delete curr->visited_nodes;
				free(curr);
			}

		} else if (is_repeat || tracker->distance > MAX_DISTANCE_FROM_MARKER) {
			// Just drop this path
			if (is_repeat) {
				printf("Track back repeat: ");
			} else {
				printf("Track back marker distance exceeded: ");
			}
			print_kmer(curr->curr_node);
			printf("\n");
			fflush(stdout);
			delete curr->visited_nodes;
			free(curr);
		} else if (curr->curr_node->fromNodes == NULL) {
			printf("Track back root: ");
			print_kmer(curr->curr_node);
			printf("\n");
			fflush(stdout);

			// We've reached a root, save it
			pthread_mutex_lock(&marker_trackback_mutex);
			roots[curr->curr_node->kmer] = curr->curr_node;
			pthread_mutex_unlock(&marker_trackback_mutex);
			delete curr->visited_nodes;
			free(curr);
		} else {

			// Add current kmer to visited kmers
			curr->visited_nodes->insert(curr->curr_node->kmer);

			// Increment distance from marker
			curr->distance += 1;

			// Move to the prev node(s)
			struct linked_node* from = curr->curr_node->fromNodes;
			char is_first = 1;

			while (from != NULL) {
				node_tracker* tracker = NULL;

				if (is_first) {
					is_first = 0;
					// Re-use the current tracker node.  Visited nodes are already tracked.
					tracker = curr;
				} else {
					// Allocate new tracker and visited nodes
					tracker = (node_tracker*) calloc(1, sizeof(node_tracker));
					tracker->visited_nodes = new sparse_hash_set<const char*, my_hash, eqstr>(*curr->visited_nodes);
					tracker->distance = curr->distance;
				}

				tracker->curr_node = from->node;
				node_stack.push(tracker);
				printf("marker traceback pushing: ");
				print_kmer(tracker->curr_node);
				printf("\n");
				from = from->next;
			}

			printf("marker traceback stack size: %d\n", node_stack.size());
			fflush(stdout);
		}
	}
}

struct marker_thread_info {
	sparse_hash_map<const char*, struct node*, my_hash, eqstr>* roots;
	struct node* node;
};

void* marker_thread(void* t) {
	marker_thread_info* info = (marker_thread_info*) t;

	get_roots_from_marker(*info->roots, info->node);

	free(info);

	pthread_mutex_lock(&running_thread_mutex);
	running_threads -= 1;
	pthread_mutex_unlock(&running_thread_mutex);
}


struct linked_node* get_roots_from_marker_nodes(struct linked_node* markers) {

	printf("Tracking back from markers\n");
	fflush(stdout);

	sparse_hash_map<const char*, struct node*, my_hash, eqstr> roots;

	pthread_t threads[MAX_THREADS];
	int num_threads = 0;
	running_threads = 0;

	// For each marker node, follow graph backwards to identify true source node(s).
	int ctr = 0;
	while (markers != NULL) {

		while (running_threads >= MAX_RUNNING_THREADS) {
			// Sleep for 10 milliseconds
			usleep(10*1000);
		}

		running_threads += 1;

		marker_thread_info* marker_info = (marker_thread_info*) calloc(1, sizeof(marker_thread_info));
		marker_info->roots = &roots;
		marker_info->node = markers->node;

		int ret = pthread_create(&(threads[num_threads]), NULL, marker_thread, marker_info);
		if (ret != 0) {
			printf("Error creating thread 1: %d\n", ret);
			fflush(stdout);
			exit(-1);
		}

		printf("Running threads: %d\n", running_threads);

		num_threads++;

//		get_roots_from_marker(roots, markers->node);
		markers = markers->next;
//		printf("Processed marker: %d\n", ctr++);
//		fflush(stdout);
	}

	while (running_threads > 0) {
		// Sleep for 10 milliseconds
		usleep(10*1000);
	}

	struct linked_node* source_nodes = NULL;

	int count = 0;

	// Iterate over all map entries to build list of source nodes
	for (sparse_hash_map<const char*, struct node*, my_hash, eqstr>::const_iterator it = roots.begin();
				 it != roots.end(); ++it) {

		const char* key = it->first;
		struct node* node = it->second;

		struct linked_node* curr;

		if (source_nodes == NULL) {
			source_nodes = (linked_node*) malloc(sizeof(linked_node));
			curr = source_nodes;
		} else {
			curr->next = (linked_node*) malloc(sizeof(linked_node));
			curr = curr->next;
		}

		curr->next = NULL;
		curr->node = node;
		count += 1;
		printf("Adding source: %d\n", count);
		fflush(stdout);
	}

	printf("Marker trace back yields: %d root nodes\n", count);
	fflush(stdout);

	return source_nodes;
}



//TODO: Order these...
struct contig {
	char* seq;
	vector<char*>* fragments;
	struct node* curr_node;
	sparse_hash_set<const char*, my_hash, eqstr>* visited_nodes;
	double score;
	int size;  // really curr_index now
	char is_repeat;
};

struct contig* new_contig() {
	struct contig* curr_contig;
	curr_contig = (contig*) calloc(1, sizeof(contig));
//	printf("seq size: %d\n", sizeof(curr_contig->seq));
//	memset(curr_contig->seq, 0, sizeof(curr_contig->seq));
	curr_contig->seq = (char*) calloc(MAX_CONTIG_SIZE, sizeof(char));
	curr_contig->size = 0;
	curr_contig->is_repeat = 0;
	curr_contig->visited_nodes = new sparse_hash_set<const char*, my_hash, eqstr>();
	curr_contig->score = 0;
	curr_contig->fragments = new vector<char*>();

	return curr_contig;
}

struct contig* copy_contig(struct contig* orig, int& num_fragments, char** all_contig_fragments) {

	// Stash any orig sequence in fragment vector so the pointer can be shared across contigs
	if (strlen(orig->seq) > 0) {
		orig->fragments->push_back(orig->seq);

		// track fragments for cleanup later
		all_contig_fragments[num_fragments++] = orig->seq;

		//TODO: Use smaller fragment block here.
		orig->seq = (char*) calloc(MAX_CONTIG_SIZE, sizeof(char));
		orig->size = 0;  // Reset index (not really size)
	}

	struct contig* copy = (contig*) calloc(sizeof(contig), sizeof(char));
	copy->seq = (char*) calloc(MAX_CONTIG_SIZE, sizeof(char));

	// Copy original fragments to new contig
	copy->fragments = new vector<char*>(*(orig->fragments));

	copy->size = orig->size;
	copy->is_repeat = orig->is_repeat;
	copy->visited_nodes = new sparse_hash_set<const char*, my_hash, eqstr>(*orig->visited_nodes);
	copy->score = orig->score;

	return copy;

}

void free_contig(struct contig* contig) {
	delete contig->visited_nodes;
	delete contig->fragments; // TODO: Indiviual fragments need to be freed!
	free(contig->seq);
	free(contig);
}

char is_node_visited(struct contig* contig, struct node* node) {
	 sparse_hash_set<const char*, my_hash, eqstr>::const_iterator it = contig->visited_nodes->find(node->kmer);
	 return it != contig->visited_nodes->end();
}

int get_contig_len(struct contig* contig) {
	int length = strlen(contig->seq);

	for (vector<char*>::iterator it = contig->fragments->begin(); it != contig->fragments->end(); ++it) {
		length += strlen(*it);
	}

	return length;
}

int output_contigs = 0;
void output_contig(struct contig* contig, int& contig_count, const char* prefix, char* contigs) {
	char buf[1024];

//	if (strlen(contigs) + strlen(contig->seq) > MAX_TOTAL_CONTIG_LEN) {
//		printf("contig string too long: %s\n", prefix);
//		exit(-1);
//	}
	pthread_mutex_lock(&contig_writer_mutex);
	output_contigs += 1;
	pthread_mutex_unlock(&contig_writer_mutex);

//	if (strlen(contig->seq) >= min_contig_length) {
	if (get_contig_len(contig) >= min_contig_length) {
		contig_count++;


		if (contig->is_repeat) {
			sprintf(buf, ">%s_%d_%f_repeat\n", prefix, output_contigs, contig->score/(double)strlen(buf));
		} else {
			sprintf(buf, ">%s_%d_%f\n", prefix, output_contigs, contig->score/(double)strlen(buf));
		}

		pthread_mutex_lock(&contig_writer_mutex);
//		fprintf(stderr, "%s%s\n", buf, contig->seq);
		fprintf(stderr, "%s", buf);
		for (vector<char*>::iterator it = contig->fragments->begin(); it != contig->fragments->end(); ++it) {
			fprintf(stderr, "%s", *it);
		}
		if (strlen(contig->seq) > 0) {
			fprintf(stderr, "%s", contig->seq);
		}
		fprintf(stderr, "\n");
		pthread_mutex_unlock(&contig_writer_mutex);

//		fflush(stderr);
	}
}




//#define OK 0
//#define TOO_MANY_PATHS_FROM_ROOT -1
//#define TOO_MANY_CONTIGS -2
//#define STOPPED_ON_REPEAT -3

#define SINK "CACACAGCCCCCAACATGCATGCTT"

int build_contigs(
		struct node* root,
		int& contig_count,
		const char* prefix,
		int max_paths_from_root,
		int max_contigs,
		char stop_on_repeat,
		char shadow_mode,
		char* contig_str) {

	int status = OK;
	stack<contig*> contigs;
//	stack<contig*> contigs_to_output;
	stack<contig*> popped_contigs;
	struct contig* root_contig = new_contig();
	root_contig->curr_node = root;
	contigs.push(root_contig);
//	vector<char*> all_contig_fragments;

	int MAX_FRAGMENTS_PER_THREAD = 10000000;
	int num_fragments = 0;
	char** all_contig_fragments = (char**) calloc(MAX_FRAGMENTS_PER_THREAD, sizeof(char*));

	int paths_from_root = 1;

	int all_contigs_len = 0;

	while ((contigs.size() > 0) && (status == OK)) {
		// Get contig from stack
		struct contig* contig = contigs.top();

//		printf ("contig: %x\n");
//		fflush(stdout);
//		printf("node: %x\n", contig->curr_node);
//		printf("kmer: %s\n", contig->curr_node->kmer);
//		fflush(stdout);
/*
		if (contig->curr_node != NULL && contig->curr_node->kmer != NULL && strncmp(contig->curr_node->kmer, SINK, kmer_size) == 0) {
			all_contigs_len += strlen(contig->seq) + 100;
			memcpy(&(contig->seq[contig->size]), contig->curr_node->kmer, kmer_size);
			output_contig(contig, contig_count, prefix, contig_str);
			free_contig(contig);
			contigs.pop();
		}
		else
*/
		if (is_node_visited(contig, contig->curr_node)) {
			printf("Repeat node: ");
			print_kmer(contig->curr_node);
			printf("\n");
			// We've encountered a repeat
			contig->is_repeat = 1;
			if ((!shadow_mode) && (!stop_on_repeat)) {
				// Add length of contig + padding for prefix
				all_contigs_len += strlen(contig->seq) + 100;
//				contigs_to_output.push(contig);

				output_contig(contig, contig_count, prefix, contig_str);
				free_contig(contig);
			} else {
				popped_contigs.push(contig);
			}
			contigs.pop();
			if (stop_on_repeat) {
				status = STOPPED_ON_REPEAT;
			}
		}
		else if (contig->curr_node->toNodes == NULL || contig->score < MIN_CONTIG_SCORE) {
//		else if (contig->curr_node->toNodes == NULL) {
			// We've reached the end of the contig.
			// Append entire current node.
			memcpy(&(contig->seq[contig->size]), contig->curr_node->kmer, kmer_size);

			// Now, write the contig
			if (!shadow_mode) {
				all_contigs_len += strlen(contig->seq) + 100;
//				contigs_to_output.push(contig);
				// only output at sink node
				output_contig(contig, contig_count, prefix, contig_str);
				free_contig(contig);
			} else {
				popped_contigs.push(contig);
			}
			contigs.pop();
//			free_contig(contig);
		}
		else {
			// Append first base from current node
			contig->seq[contig->size++] = contig->curr_node->kmer[0];
			if (contig->size >= MAX_CONTIG_SIZE) {
				char kmer[1024];
				memset(kmer, 0, 1024);
				strncpy(kmer, contig->curr_node->kmer, kmer_size);
				printf("Max contig size exceeded at node: %s\n", kmer);

				//TODO: Provide different status
				status = TOO_MANY_CONTIGS;
				break;
			}

			contig->visited_nodes->insert(contig->curr_node->kmer);

			//			struct linked_node* selected_to_node = NULL;
			//			int max_freq = -1;
			//
			//			while (to_linked_node != NULL) {
			//
			//				struct node* to_node = to_linked_node->node;
			//				if (to_node->frequency > max_freq) {
			//					max_freq = to_node->frequency;
			//					selected_to_node = to_linked_node;
			//				}
			//
			//				to_linked_node = to_linked_node->next;
			//			}
			//
			//			contig->curr_node = selected_to_node->node;
			//			paths_from_root++;


			// Count total edges
			int total_edge_count = 0;
			struct linked_node* to = contig->curr_node->toNodes;

			while (to != NULL) {
				total_edge_count = total_edge_count + to->node->frequency;
				to = to->next;
			}

			double log10_total_edge_count = log10(total_edge_count);

			// Move current contig to next "to" node.
			struct linked_node* to_linked_node = contig->curr_node->toNodes;

			contig->curr_node = to_linked_node->node;
			paths_from_root++;

			// If there are multiple "to" nodes, branch the contig and push on stack
			to_linked_node = to_linked_node->next;

			while (to_linked_node != NULL) {
				struct contig* contig_branch = copy_contig(contig, num_fragments, all_contig_fragments);
				contig_branch->curr_node = to_linked_node->node;
				contig_branch->score = contig_branch->score + log10(contig_branch->curr_node->frequency) - log10_total_edge_count;
				contigs.push(contig_branch);
				to_linked_node = to_linked_node->next;
				paths_from_root++;
			}

			contig->score = contig->score + log10(contig->curr_node->frequency) - log10(total_edge_count);
		}

		if (contig_count >= max_contigs) {
			status = TOO_MANY_CONTIGS;
		}

		if (paths_from_root >= max_paths_from_root) {
			status = TOO_MANY_PATHS_FROM_ROOT;
		}
	}
/*
	if (status == OK) {

		while (contigs_to_output.size() > 0) {
			struct contig* contig = contigs_to_output.top();
			output_contig(contig, contig_count, prefix, contig_str);

			contigs_to_output.pop();
			free_contig(contig);
		}
	}


	// Cleanup stranded contigs in case processing stopped.
	while (contigs_to_output.size() > 0) {
		struct contig* contig = contigs_to_output.top();

		contigs_to_output.pop();
		free_contig(contig);
	}
*/
	while (contigs.size() > 0) {
		struct contig* contig = contigs.top();
		contigs.pop();
		free_contig(contig);
	}

	while (popped_contigs.size() > 0) {
		struct contig* contig = popped_contigs.top();
		popped_contigs.pop();
		free_contig(contig);
	}

	printf("Freeing %d fragments\n", num_fragments);
	for (int i=0; i<num_fragments; i++) {
		free(all_contig_fragments[i]);
	}

	free(all_contig_fragments);

	// Cleanup contig fragments
//	for (vector<char*>::iterator it = all_contig_fragments.begin(); it != all_contig_fragments.end(); ++it) {
//		free(*it);
//	}

	return status;
}

//TODO: Use thread pool instead of spawning threads
void* worker_thread(void* t) {

	struct linked_node* root = (struct linked_node*) t;

	while (root != NULL) {

		struct node* source = root->node;

		int contig_count;
		const char* prefix = "foo";
		int max_paths_from_root = 500000000;
		int max_contigs = 50000000;
		char stop_on_repeat = false;
		char shadow_mode = false;
		char* contig_str = NULL;
		int status = build_contigs(source, contig_count, prefix, max_paths_from_root, max_contigs, stop_on_repeat,
				shadow_mode, contig_str);

		if (status != OK) {
			printf("Status: %d\n", status);
			fflush(stdout);
			exit(status);
		}

		root = root->next;
	}

	pthread_mutex_lock(&running_thread_mutex);
	running_threads -= 1;
	pthread_mutex_unlock(&running_thread_mutex);
}

struct simple_node {
	char* seq;
	vector<char*>* to_kmers;
	int id;
	unsigned short freq1;
	unsigned short freq2;
};

int is_simplify_start_point(struct node* node) {
	int is_start = 0;

	// Root nodes and nodes with multiple incoming edges are
	// eligible starting points for simplification
	if (node->fromNodes == NULL) {
		is_start = 1;
	} else {
		linked_node* from = node->fromNodes;
		if (from->next != NULL) {
			is_start = 1;
		}
	}

	// Nodes with a from node containing multiple outgoing edges
	// are a simplification start point
	if (!is_start) {
		linked_node* from = node->fromNodes;

		if (from->node->toNodes != NULL && from->node->toNodes->next != NULL) {
			is_start = 1;
		}
	}

	return is_start;
}

int is_simplify_stop_point(struct node* node) {
	int is_stop = 0;

	// We've reached a stop point if the number of
	// outgoing edges != 1
	if (node->toNodes == NULL) {
		is_stop = 1;
	} else {
		linked_node* to = node->toNodes;
		if (to->next != NULL) {
			is_stop = 1;
		}
	}

	// If the # of outgoing edges is 1 and the to node
	// has multiple incoming edges, this is a stop point
	if (!is_stop) {
		linked_node* to = node->toNodes;

		if (to->node->fromNodes != NULL && to->node->fromNodes->next != NULL) {
			is_stop = 1;
		}
	}

	return is_stop;
}

sparse_hash_map<const char*, struct simple_node*, my_hash, eqstr>* simplify_graph(sparse_hash_map<const char*, struct node*, my_hash, eqstr>* nodes) {

	sparse_hash_map<const char*, struct simple_node*, my_hash, eqstr>* simple_nodes = new sparse_hash_map<const char*, struct simple_node*, my_hash, eqstr>();

	int id = 0;
	int is_in_linear_path = 0;
	vector<struct node*> curr_nodes;
	for (sparse_hash_map<const char*, struct node*, my_hash, eqstr>::const_iterator it = nodes->begin();
				 it != nodes->end(); ++it) {

		const char* key = it->first;
		struct node* node = it->second;

		if (is_simplify_start_point(node)) {
			printf("start point:\n");
			print_node(node);
			printf("-----\n");
			is_in_linear_path = 1;
			curr_nodes.clear();
			curr_nodes.push_back(node);

			struct linked_node* to_node = node->toNodes;
			// Iterate over all nodes in linear path (skipping singletons)
			while (to_node != NULL) {
				node = to_node->node;

//				printf("path:\n");
//				print_node(node);
//				printf("\n----\n");

				if (is_simplify_stop_point(node)) {
					int len = curr_nodes.size();
					if (node->toNodes == NULL) {
						len += kmer_size;
					} else {
						len += 1;
					}

					char* seq = (char*) calloc(len+1, sizeof(char));
					simple_node* compressed_node = (simple_node*) calloc(1, sizeof(simple_node));
					compressed_node->seq = seq;
					compressed_node->id = id++;
					compressed_node->to_kmers = new vector<char*>();

	//				printf("here\n");
	//				fflush(stdout);
	//				printf("to_kmers size: %d\n", compressed_node->to_kmers.size());
	//				fflush(stdout);

					int idx = 0;

					// Extract sequence from nodes in linear path
					for (vector<struct node*>::iterator node_it = curr_nodes.begin(); node_it != curr_nodes.end(); node_it++) {
						seq[idx++] = (*node_it)->kmer[0];
					}

					printf("to nodes: %x\n", node->toNodes);
					if (node->toNodes == NULL) {
						// Last node in path.  Append entire kmer
						memcpy(&seq[idx], node->kmer, kmer_size);
					} else {
						// We've reached a fork in the graph.  Record the to_kmers.
						seq[idx] = node->kmer[0];
						struct linked_node* to_node =  node->toNodes;
						while (to_node != NULL) {
							printf("pushing...\n");
							compressed_node->to_kmers->push_back(to_node->node->kmer);
							to_node = to_node->next;
						}
					}

					curr_nodes.push_back(node);

					compressed_node->freq1 = curr_nodes.front()->frequency;
					compressed_node->freq2 = node->frequency;

					printf("seq: %s\n", compressed_node->seq);

					// Use the first kmer contributing to this compressed node as the key into the storage map.
					if (curr_nodes.front()->fromNodes == NULL && node->toNodes == NULL) {
						// This is an isolated path.  Keep only if longer than read length.
						if (strlen(compressed_node->seq) > read_length) {
							(*simple_nodes)[curr_nodes.front()->kmer] = compressed_node;
						}
					} else {
						(*simple_nodes)[curr_nodes.front()->kmer] = compressed_node;
					}

					is_in_linear_path = 0;
					curr_nodes.clear();

					to_node = NULL;

				} else {
					curr_nodes.push_back(node);
					// There must be a single to node at this point.
					// Advance to next node.
					to_node = node->toNodes;
				}
			}
		}
	}

	return simple_nodes;
}

void cleanup_simple_graph(sparse_hash_map<const char*, struct simple_node*, my_hash, eqstr>* simple_nodes) {
	for (sparse_hash_map<const char*, struct simple_node*, my_hash, eqstr>::const_iterator it = simple_nodes->begin();
				 it != simple_nodes->end(); ++it) {

		const char* kmer = it->first;
		struct simple_node* curr_node = it->second;

//		simple_nodes->erase(kmer);
		free(curr_node->seq);
		delete curr_node->to_kmers;
		free(curr_node);
	}

	free(simple_nodes);
}

void get_seq_str(char* seq, char* buf, int buf_size) {
	char temp[128];

	memset(buf, 0, buf_size);
	memset(temp, 0, 128);

	int seq_len = strlen(seq);
	if (seq_len > 25) {
		for (int i=0; i<10; i++) {
			temp[i] = seq[i];
			temp[13+i] = seq[seq_len-(10-i)-1];
		}
		temp[10] = '.';
		temp[11] = '.';
		temp[12] = '.';
	} else {
		strcpy(temp, seq);
	}

	sprintf(buf, "%s (%d)", temp, strlen(seq));
}

void write_graph(sparse_hash_map<const char*, struct node*, my_hash, eqstr>* nodes) {

	sparse_hash_map<const char*, struct simple_node*, my_hash, eqstr>* simple_nodes = simplify_graph(nodes);

	const char* filename = "graph.dot";

	FILE *fp = fopen(filename, "w");

	fprintf(fp, "digraph assemblyGraph {\n");

	for (sparse_hash_map<const char*, struct simple_node*, my_hash, eqstr>::const_iterator it = simple_nodes->begin();
				 it != simple_nodes->end(); ++it) {
		const char* kmer = it->first;
		struct simple_node* curr_node = it->second;

		printf("node: %x\n", curr_node);
		fflush(stdout);
		printf("to_kmers: %x\n", curr_node->to_kmers);
		fflush(stdout);
		printf("to_kmers: %d\n", curr_node->to_kmers->size());
		fflush(stdout);

		// Output vertex
		char sequence[128];
		get_seq_str(curr_node->seq, sequence, 128);
		fprintf(fp, "\tvertex_%d [label=\"%s\",shape=box]\n", curr_node->id, sequence);

		// Output outgoing edges
		for (vector<char*>::iterator kmer_it = curr_node->to_kmers->begin(); kmer_it != curr_node->to_kmers->end(); kmer_it++) {
			char* kmer = (*kmer_it);

			printf("from_seq: %s\n", curr_node->seq);
			fflush(stdout);
			printf("to_kmer: \n");
			print_kmer(kmer);

			print_node((*nodes)[kmer]);

			printf("\n");
			fflush(stdout);

			if (simple_nodes->find(kmer) != simple_nodes->end()) {

				struct simple_node* to_node = (*simple_nodes)[kmer];
				printf("edge to: %x\n", to_node);
				fflush(stdout);


				fprintf(fp, "\tvertex_%d -> vertex_%d [label=\"%d,%d\"];", curr_node->id, to_node->id,
											curr_node->freq2, to_node->freq1);
				/*
				if (to_node != NULL) {
					fprintf(fp, "\tvertex_%d -> vertex_%d [label=\"%d,%d\"];", curr_node->id, to_node->id,
							curr_node->freq2, to_node->freq1);
				} else {
					printf("Missing to node: \n");
					print_kmer(kmer);
					printf("\n-----------\n");
					fflush(stdout);
				}
				*/
			} else {
				printf("Missing to node: ");
				print_kmer(kmer);
				printf("\n-----------\n");
				fflush(stdout);
			}
		}
	}

	fprintf(fp,"}");

	cleanup_simple_graph(simple_nodes);
}

void cleanup(sparse_hash_map<const char*, struct node*, my_hash, eqstr>* nodes, struct struct_pool* pool) {

	// Free linked lists
	for (sparse_hash_map<const char*, struct node*, my_hash, eqstr>::const_iterator it = nodes->begin();
	         it != nodes->end(); ++it) {
		struct node* node = it->second;

		if (node != NULL) {
			cleanup(node->toNodes);
			cleanup(node->fromNodes);
		}
	}

	for (int i=0; i<=pool->node_pool->block_idx; i++) {
		free(pool->node_pool->nodes[i]);
	}

	free(pool->node_pool->nodes);
	free(pool->node_pool);

	for (int i=0; i<=pool->read_pool->block_idx; i++) {
		free(pool->read_pool->reads[i]);
	}

	free(pool->read_pool->reads);
	free(pool->read_pool);

	free(pool);
}

void dump_graph(sparse_hash_map<const char*, struct node*, my_hash, eqstr>* nodes, const char* filename) {

	FILE* fp = fopen(filename, "w");

	for (sparse_hash_map<const char*, struct node*, my_hash, eqstr>::const_iterator it = nodes->begin();
				 it != nodes->end(); ++it) {

		const char* key = it->first;
		node* curr_node = it->second;

		char output[1024];
		memset(output, 0, 1024);
		strncpy(output, curr_node->kmer, kmer_size);
		int idx = kmer_size;
		strncpy(&(output[idx]), "->", 2);
		idx += 2;

		if (curr_node != NULL) {
			////////////////////////////////////////////////
			// Check to node list for low frequency edges
			struct linked_node* to_node = curr_node->toNodes;

			while (to_node != NULL) {
				strncpy(&(output[idx]), to_node->node->kmer, kmer_size);
				idx += kmer_size;
				output[idx] = ',';
				idx += 1;
				to_node = to_node->next;
			}
		}

		fprintf(fp, "%s\n", output);
	}

	fclose(fp);
}


char* assemble(const char* input,
			   const char* unaligned_input,
			  const char* output,
			  const char* prefix,
			  int truncate_on_repeat,
			  int max_contigs,
			  int max_paths_from_root,
			  int input_read_length,
			  int input_kmer_size,
			  const char* bcr_fasta) {


	read_length = input_read_length;

//	min_contig_length = read_length + 1;

	//TODO: Parameterize mcl - shorter for unaligned region?
/*
	if (truncate_on_repeat) {
		min_contig_length = read_length + 1;
	} else {
		min_contig_length = 150;
	}
*/

	kmer_size = input_kmer_size;

	printf("Init root scoring start\n");
	init_root_scoring(bcr_fasta, kmer_size);
	printf("Init root scoring done\n");

	struct struct_pool* pool = init_pool();
	sparse_hash_map<const char*, struct node*, my_hash, eqstr>* nodes = new sparse_hash_map<const char*, struct node*, my_hash, eqstr>();

	long startTime = time(NULL);
	printf("Assembling: -> %s\n", output);
	nodes->set_deleted_key(NULL);

	pthread_mutex_init(&running_thread_mutex, NULL);
	pthread_mutex_init(&contig_writer_mutex, NULL);
	pthread_mutex_init(&marker_trackback_mutex, NULL);

	// Init scoring matrix
	init_score_seq(10000,200);

	build_graph2(input, nodes, pool, 1);

	char isUnalignedRegion = !truncate_on_repeat;
	printf("Prune graph 1...\n");
	fflush(stdout);
	prune_graph(nodes, isUnalignedRegion);

	struct linked_node* root_nodes = NULL;
	root_nodes = identify_root_nodes(nodes);

	if (unaligned_input != NULL) {
		build_graph2(unaligned_input, nodes, pool, 0);
	}

	int status = -1;

	if (nodes->size() >= MAX_NODES) {
		status = TOO_MANY_NODES;
		printf("Graph too complex for region: %s\n", prefix);
	}

	//TODO: Set this explicitly
//	char isUnalignedRegion = !truncate_on_repeat;
	printf("Prune graph 2...\n");
	fflush(stdout);
	prune_graph(nodes, isUnalignedRegion);
	printf("Prune graph 2 Done...\n");
	fflush(stdout);

	dump_graph(nodes, "graph.txt");


	// The initial set of root nodes are used as markers (or breadcrumbs)
	// We start from a marker and follow the trail back to the true source nodes
//	root_nodes = get_roots_from_marker_nodes(root_nodes);

//	if (status != TOO_MANY_NODES) {
//		root_nodes = identify_root_nodes(nodes);
//	}

	int contig_count = 0;
	char truncate_output = 0;

	char* contig_str = (char*) malloc(MAX_TOTAL_CONTIG_LEN);
	memset(contig_str, 0, MAX_TOTAL_CONTIG_LEN);


	pthread_t threads[MAX_THREADS];
	int num_threads = 0;
	running_threads = 0;

//	FILE *fp = fopen(output, "w");
	int num_roots = 0;
	pthread_mutex_init(&running_thread_mutex, NULL);
	pthread_mutex_init(&contig_writer_mutex, NULL);

	linked_node* root_link = NULL;

	while (root_nodes != NULL) {

//		int shadow_count = 0;

		while (running_threads >= MAX_RUNNING_THREADS) {
			// Sleep for 10 milliseconds
			usleep(10*1000);
		}

		int ROOT_NODES_PER_THREADS = 10;
//		int ROOT_NODES_PER_THREADS = 2;

		int root_count = 1;
		linked_node* next_root_start = root_nodes;
		linked_node* last = root_nodes;
		while (root_nodes != NULL && root_count < ROOT_NODES_PER_THREADS) {
			root_nodes = root_nodes->next;
			if (root_nodes != NULL) {
				last = root_nodes;
			}
			num_roots++;
			root_count += 1;
		}

		root_nodes = last->next;
		last->next = NULL;

		running_threads += 1;
//		int ret = pthread_create(&(threads[num_threads]), NULL, worker_thread, root_nodes->node);
		int ret = pthread_create(&(threads[num_threads]), NULL, worker_thread, next_root_start);


		if (ret != 0) {
			printf("Error creating thread 1: %d\n", ret);
			fflush(stdout);
			exit(-1);
		}

		printf("Running threads: %d\n", running_threads);

		num_threads++;

/*
		status = build_contigs(root_nodes->node, contig_count, prefix, max_paths_from_root, max_contigs, truncate_on_repeat, false, contig_str);

		switch(status) {
			case TOO_MANY_CONTIGS:
				printf("TOO_MANY_CONTIGS: %s\n", prefix);
				contig_count = 0;
				break;
			case STOPPED_ON_REPEAT:
				printf("STOPPED_ON_REPEAT: %s\n", prefix);
				contig_count = 0;
				break;
			case TOO_MANY_PATHS_FROM_ROOT:
				char kmer[1024];
				memset(kmer, 0, 1024);
				strncpy(kmer, root_nodes->node->kmer, kmer_size);
				printf("TOO_MANY_PATHS_FROM_ROOT: %s - %s\n", prefix, kmer);
				break;
		}

		// If too many contigs or abort due to repeat, break out of loop and truncate output.
		if ((status == TOO_MANY_CONTIGS) || (status == STOPPED_ON_REPEAT)) {
			truncate_output = 1;
			break;
		}
*/

//S		root_nodes = root_nodes->next;

		if ((num_roots % 10) == 0) {
			printf("Processed %d root nodes\n", num_roots);
			fflush(stdout);
		}
	}

	// Wait for all threads to complete
	while (running_threads > 0) {
		// Sleep for 10 milliseconds
		usleep(10*1000);
	}

	pthread_mutex_destroy(&running_thread_mutex);
	pthread_mutex_destroy(&contig_writer_mutex);
	pthread_mutex_destroy(&marker_trackback_mutex);

	printf("output contigs: %d\n", output_contigs);

//	write_graph(nodes);

	cleanup(nodes, pool);

	delete nodes;

	long stopTime = time(NULL);

	if (kmer_size != input_kmer_size) {
		printf("What!!?? %d : %d\n", kmer_size, input_kmer_size);
	}
	assert(kmer_size == input_kmer_size);
	printf("Done assembling(%ld): %s, %d\n", (stopTime-startTime), output, contig_count);

	if (status == OK || status == TOO_MANY_PATHS_FROM_ROOT) {
		return contig_str;
	} else if (status == STOPPED_ON_REPEAT) {
		strcpy(contig_str, "<REPEAT>");
		return contig_str;
	} else {
		strcpy(contig_str, "<ERROR>");
		return contig_str;
	}
}

/*
extern "C"
 JNIEXPORT jstring JNICALL Java_abra_NativeAssembler_assemble
   (JNIEnv *env, jobject obj, jstring j_input, jstring j_output, jstring j_prefix,
    jint j_truncate_on_output, jint j_max_contigs, jint j_max_paths_from_root,
    jint j_read_length, jint j_kmer_size, jint j_min_node_freq, jint j_min_base_quality)
 {

     //Get the native string from javaString
     //const char *nativeString = env->GetStringUTFChars(javaString, 0);
	const char* input  = env->GetStringUTFChars(j_input, 0);
	const char* output = env->GetStringUTFChars(j_output, 0);
	const char* prefix = env->GetStringUTFChars(j_prefix, 0);
	int truncate_on_output = j_truncate_on_output;
	int max_contigs = j_max_contigs;
	int max_paths_from_root = j_max_paths_from_root;
	int read_length = j_read_length;
	int kmer_size = j_kmer_size;
	min_node_freq = j_min_node_freq;
	min_base_quality = j_min_base_quality;

	printf("Abra JNI entry point v0.79\n");

	printf("input len: %s : %d\n", prefix, strlen(input));
	printf("output: %s\n", output);
	printf("prefix: %s\n", prefix);
	printf("truncate_on_output: %d\n", truncate_on_output);
	printf("max_contigs: %d\n", max_contigs);
	printf("max_paths_from_root: %d\n", max_paths_from_root);
	printf("read_length: %d\n", read_length);
	printf("kmer_size: %d\n", kmer_size);
	printf("min node freq: %d\n", min_node_freq);
	printf("min base quality: %d\n", min_base_quality);

	char* contig_str = assemble(input, NULL, output, prefix, truncate_on_output, max_contigs, max_paths_from_root, read_length, kmer_size, "bcr_fasta here");
	jstring ret = env->NewStringUTF(contig_str);

     //DON'T FORGET THIS LINE!!!
    env->ReleaseStringUTFChars(j_input, input);
    env->ReleaseStringUTFChars(j_output, output);
    env->ReleaseStringUTFChars(j_prefix, prefix);
    free(contig_str);

    fflush(stdout);

    return ret;
 }
 */


char* load_file(const char* filename) {
	FILE    *infile;
	char    *buffer;
	long    numbytes;

	printf("Loading: %s\n", filename);
	fflush(stdout);

	infile = fopen(filename, "r");
	fseek(infile, 0L, SEEK_END);
	numbytes = ftell(infile);
	printf("numbytes: %ld\n", numbytes);
	fflush(stdout);
	buffer = (char*) malloc(numbytes+1);
	memset(buffer, 0, numbytes+1);
	fseek(infile, 0L, SEEK_SET);
	fread(buffer, sizeof(char), numbytes, infile);
	fclose(infile);
	return buffer;
}

int main(int argc, char* argv[]) {

	char* input_file = argv[1];
	char* unaligned_input_file = NULL;

    min_node_freq = 5;
    min_base_quality = 150;

	//unaligned_input_file = argv[2];
	min_node_freq = atoi(argv[2]);
	min_base_quality = atoi(argv[3]);
//		min_base_quality = min_node_freq * 30;
	if (min_base_quality >= MAX_QUAL_SUM) {
		min_base_quality = MAX_QUAL_SUM-1;
	}

	int read_len = atoi(argv[4]);

	MIN_CONTIG_SCORE = atof(argv[5]);

	MAX_RUNNING_THREADS = atoi(argv[6]);

	min_contig_length = 501;

//	char* input = load_file("/datastore/nextgenout4/seqware-analysis/lmose/vdj/TCGA-FF-8041-01A-11R-2213-07/igkv4_1.b.reads");
//	char* unaligned_input = load_file("/datastore/nextgenout4/seqware-analysis/lmose/vdj/TCGA-FF-8041-01A-11R-2213-07/igkv4_1_unaligned.c.reads");

//	char* input = load_file("/datastore/nextgenout4/seqware-analysis/lmose/vdj/TCGA-FF-8041-01A-11R-2213-07/pass2/pass2.reads");
//	char* unaligned_input = load_file("/datastore/nextgenout4/seqware-analysis/lmose/vdj/TCGA-FF-8041-01A-11R-2213-07/pass2/pass2_unaligned.reads");

	//char* input = load_file("/datastore/nextgenout4/seqware-analysis/lmose/vdj/TCGA-FF-8041-01A-11R-2213-07/target.reads");
//	char* input = load_file("/home/lmose/dev/vdj/viz/test1.reads");
//	char* input = load_file("/home/lmose/dev/vdj/viz/test2.reads");

//	char* input = load_file("/datastore/nextgenout4/seqware-analysis/lmose/vdj/mouse/test1/pass2/rna.vdj.reads");
	char* input = load_file(input_file);
	char* unaligned_input = NULL;

	if (unaligned_input_file != NULL) {
		unaligned_input = load_file(unaligned_input_file);
	}

	char* bcr_fasta = "/datastore/nextgenout4/seqware-analysis/lmose/vdj/mouse/mm10.bcr.constant.fa";

//	char* unaligned_input = load_file("/datastore/nextgenout4/seqware-analysis/lmose/vdj/TCGA-FF-8041-01A-11R-2213-07/all_unaligned.reads");

        char* output = assemble(
//                "/datastore/nextgenout4/seqware-analysis/lmose/platinum/test/mtest.reads",
//                "/datastore/nextgenout4/seqware-analysis/lmose/platinum/test/mtest.fa",
//                "/datastore/nextgenout4/seqware-analysis/lmose/platinum/test/ftest.reads",
//                "/datastore/nextgenout4/seqware-analysis/lmose/platinum/test/ftest.fa",
//		"/datastore/nextgenout4/seqware-analysis/lmose/CALGB_40603/data/CALGB03-1049119-I-31/unaligned.reads",
		input,
		unaligned_input,
		"",
                "foo",
                false,
                50000000,
                500000000,
                read_len,
                25,
                bcr_fasta);

}


