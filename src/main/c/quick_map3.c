#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <algorithm>
#include <utility>
#include <sparsehash/sparse_hash_map>

#include "hash_utils.h"
#include "quick_map3.h"

using namespace std;

using google::sparse_hash_map;

//int READ_LEN = 50;
//int MIN_INSERT = 180 - 60; // 120
//int MAX_INSERT = 180 + 60; // 240

// Global from assembler2_vdj.c
extern int read_length;
int READ_LEN = 0;
int MIN_INSERT = 50;
int MAX_INSERT = 400;

//
// String comparison bounded at READ_LEN
struct qm_eqstr
{
  bool operator()(const char* s1, const char* s2) const
  {
    return (s1 == s2) || (s1 && s2 && strncmp(s1, s2, READ_LEN) == 0);
  }
};

//
// String hash bounded at READ_LEN
struct qm_hash
{
	uint64_t operator()(const char* seq) const
	{
		return MurmurHash64A(seq, READ_LEN, 97);
	}
};

// Allocate 1GB at a time
#define READ_BLOCK 1000000000

char* read_buf;
char* read_buf_start;
char* qual_buf;
char* qual_buf_start;

#define MAX_LINE 1024
#define MAX_CONTIG_LEN 10000

//
// Key = read sequence, Value = vector of read_info
sparse_hash_map<const char*, struct read_vec*, qm_hash, qm_eqstr>* reads = new sparse_hash_map<const char*, struct read_vec*, qm_hash, qm_eqstr>();

void quick_map_init() {
	READ_LEN = read_length;
	reads->set_deleted_key(NULL);
}

void advance_read_buf() {
	// Advance read buffer to next open slot (TODO: Better to stay on word boundary?)
	read_buf = read_buf + strlen(read_buf) + 1;

	// If we are approaching the end of the read buffer, allocate more space
	if (read_buf-read_buf_start > READ_BLOCK-1024) {
		fprintf(stderr, "increasing read buf\n");
		fflush(stderr);
		read_buf = (char*) calloc(READ_BLOCK, sizeof(char));
		read_buf_start = read_buf;
	}
}

//TODO: Factor out common code
void advance_qual_buf() {
	// Advance read buffer to next open slot (TODO: Better to stay on word boundary?)
	qual_buf = qual_buf + strlen(qual_buf) + 1;

	// If we are approaching the end of the read buffer, allocate more space
	if (qual_buf-qual_buf_start > READ_BLOCK-1024) {
		fprintf(stderr, "increasing qual buf\n");
		fflush(stderr);
		qual_buf = (char*) calloc(READ_BLOCK, sizeof(char));
		qual_buf_start = qual_buf;
	}
}

/*
char complement(char ch) {
	switch(ch) {
		case 'A':
			return 'T';
		case 'T':
			return 'A';
		case 'C':
			return 'G';
		case 'G':
			return 'C';
		default:
			return ch;
	}
}

int rc(char* input, char* output) {
	int out_idx = 0;
	for (int i=strlen(input)-1; i >=0; i--) {
		output[out_idx++] = complement(input[i]);
	}
	output[strlen(input)] = '\0';
}

int reverse(char* input, char* output) {
	int out_idx = 0;
	for (int i=strlen(input)-1; i >=0; i--) {
		output[out_idx++] = input[i];
	}
	output[strlen(input)] = '\0';
}
*/

void add_read_info(char* read_id, char* seq, char* quals, char read_num, char is_rc) {

	read_vec* seq_reads = (*reads)[seq];

	// No hit in hash map.  Add new entry
	if (seq_reads == NULL) {
		seq_reads = (read_vec*) calloc(1, sizeof(read_vec));
		seq_reads->reads = new vector<read_info*>();
		seq_reads->seq = seq;
		(*reads)[seq] = seq_reads;
	}

	read_info* read_info1 = (read_info*) calloc(1, sizeof(read_info));

	read_info1->id = read_id;
	read_info1->read_num = read_num;
	read_info1->is_rc = is_rc;
	read_info1->seq = seq_reads->seq;
	read_info1->quals = quals;

	seq_reads->reads->push_back(read_info1);

//	printf("seq: %s , read_info1->read_num: %d\n", seq, read_info1->read_num);
}

//TODO: Output base qualities
void output_mapping(char* contig_id, map_info* r1, map_info* r2, int insert) {

	char* read_id;

	if (r1->info->id[0] == '@') {
		read_id = &(r1->info->id[1]);
	} else {
		read_id = r1->info->id;
	}

	int flag1 = 0x1l | 0x2 |  (r1->info->is_rc ? 0x10 : 0x20) | 0x40;
	int flag2 = 0x1l | 0x2 |  (r2->info->is_rc ? 0x10 : 0x20) | 0x80;

	char seq[256];
	char quals[256];

	const char* format = "%s\t%d\t%s\t%d\t255\t%dM\t=\t%d\t%d\t%s\t%s\n";

	// Only need to add null terminators once.
	seq[READ_LEN] = '\0';
	quals[READ_LEN] = '\0';

	strncpy(seq, r1->info->seq, READ_LEN);
	strncpy(quals, r1->info->quals, READ_LEN);
	printf(format, read_id, flag1, contig_id, r1->pos, READ_LEN, r2->pos, insert, seq, quals);

	strncpy(seq, r2->info->seq, READ_LEN);
	strncpy(quals, r2->info->quals, READ_LEN);
	printf(format, read_id, flag2, contig_id, r2->pos, READ_LEN, r1->pos, insert, seq, quals);
}

char contains_read(char* read) {
	sparse_hash_map<const char*, struct read_vec*, qm_hash, qm_eqstr>::const_iterator it = reads->find(read);
	return it != reads->end();
}

void quick_map_process_contig(char* contig_id, char* contig, vector<mapped_pair>& mapped_reads,
		vector<pair<int, int> >& start_positions, char should_output) {

	int read1_count = 0;
	// TODO: Grow dynamically
	int MAX_READ_PAIRS = 10000000;
	map_info** read1 = (map_info**) calloc(MAX_READ_PAIRS, sizeof(map_info*));

	sparse_hash_map<const char*, struct map_info*, vjf_hash, vjf_eqstr> read2;

	// Load read 1 matches into vector
	// Load read 2 matches into map
	for (int i=0; i<strlen(contig)-READ_LEN; i++) {
		if (contains_read(contig+i)) {
			read_vec* read_v = (*reads)[contig+i];
			for (vector<read_info*>::iterator it = read_v->reads->begin(); it != read_v->reads->end(); ++it) {
				read_info* r_info = *it;

				map_info* m_info = (map_info*) calloc(1, sizeof(map_info));
				m_info->info = r_info;
				m_info->pos = i + 1;

				if (r_info->read_num == 1) {
					read1[read1_count++] = m_info;
				} else {
					// TODO: Handle read 2 multi-mappers
					read2[r_info->id] = m_info;
				}
			}
		}
	}

	// Go through read 1 array looking for matches in read 2 map
	for (int i=0; i<read1_count; i++) {

		map_info* r1 = read1[i];
		map_info* r2 = read2[r1->info->id];

		if (r2 != NULL) {

			if (r1->info->is_rc != r2->info->is_rc) {
				short insert = abs(r1->pos - r2->pos) + READ_LEN;
				if (insert >= MIN_INSERT && insert <= MAX_INSERT) {
					// We have a hit.  Output
					mapped_pair read_pair;
					read_pair.contig_id = contig_id;
					read_pair.r1 = r1->info;
					read_pair.r2 = r2->info;
					read_pair.pos1 = r1->pos;
					read_pair.pos2 = r2->pos;
					read_pair.insert = insert;
					mapped_reads.push_back(read_pair);
					start_positions.push_back(make_pair((int) read_pair.pos1, (int) read_pair.pos2));
					start_positions.push_back(make_pair((int) read_pair.pos2, (int) read_pair.pos1));
					if (should_output) {
						output_mapping(contig_id, r1, r2, insert);
					}
				}
			}
		}
	}

	// Now sort the start positions
	sort(start_positions.begin(), start_positions.end());

	// Clean up memory
	for (int i=0; i<read1_count; i++) {
		free(read1[i]);
	}

	free(read1);

	for (sparse_hash_map<const char*, struct map_info*, vjf_hash, vjf_eqstr>::const_iterator it = read2.begin();
			it != read2.end(); ++it) {

		map_info* m_info = it->second;
		free(m_info);
	}
}

void quick_map_process_contig(char* contig_id, char* contig, vector<mapped_pair>& mapped_reads,
		vector<pair<int, int> >& start_positions) {

	quick_map_process_contig(contig_id, contig, mapped_reads, start_positions, 0);
}

void output_header(char* contig_file) {
	fprintf(stderr, "Writing header\n");
	FILE* fp = fopen(contig_file, "r");
	char contig_id[256];
	char* contig = (char*) calloc(MAX_CONTIG_LEN, sizeof(char));

	printf("@HD\tVN:1.4\tSO:unsorted\n");

	while (fgets(contig, MAX_CONTIG_LEN, fp) != NULL) {
		contig[strlen(contig)-1] = '\0';  // strip newline

		if (contig[0] == '>') {
			strncpy(contig_id, &(contig[1]), 256);
		} else {
			printf("@SQ\tSN:%s\tLN:%d\n", contig_id, strlen(contig));
		}
	}

	//TODO: Add @PG
	/*
	char pg[10000];
	pg[0] = '\0';
	int chars_left = 9999;
	for (int i=0; i<argc; i++) {
		strncat(pg, argv[i], chars_left);
		strcat(pg, " ");
		chars_left = 10000 - strlen(pg);
	}

	printf("@PG\tID:quickmap\tPN:quickmap\tCL:%s\n", pg);
	*/

	fclose(fp);
	fprintf(stderr, "Done writing header\n");
	fflush(stderr);
}

void quick_map_process_contig_file(char* contig_file) {
	output_header(contig_file);

	FILE* fp = fopen(contig_file, "r");
	char contig_id[256];
	char* contig = (char*) calloc(MAX_CONTIG_LEN, sizeof(char));
	int num_contigs = 0;

	while (fgets(contig, MAX_CONTIG_LEN, fp) != NULL) {
		contig[strlen(contig)-1] = '\0';

		if (contig[0] == '>') {
			strncpy(contig_id, &(contig[1]), 256);
		} else {
//			fprintf(stderr, "Processing contig [%s]\n", contig_id);
//			fflush(stderr);

			vector<mapped_pair> mapped_reads;
			vector<pair<int,int> > start_positions;
			quick_map_process_contig(contig_id, contig, mapped_reads, start_positions, 1);
			num_contigs++;
			if ((num_contigs % 1000) == 0) {
				fprintf(stderr, "[%d] contigs processed\n", num_contigs);
				fflush(stderr);
				fflush(stdout);
			}
		}
	}
	fclose(fp);
}

/*
int main(int argc, char** argv) {
	char* fastq1 = argv[1];
	char* fastq2 = argv[2];
	char* contigs = argv[3];
	READ_LEN = atoi(argv[4]);
	MIN_INSERT = atoi(argv[5]);
	MAX_INSERT = atoi(argv[6]);

	fprintf(stderr, "fastq1: %s\n", fastq1);
	fprintf(stderr, "fastq2: %s\n", fastq2);
	fprintf(stderr, "contigs: %s\n", contigs);

	load_reads(fastq1, fastq2);
	output_header(contigs, argc, argv);
	process_contigs(contigs);

//	for (sparse_hash_map<const char*, struct read_vec*, my_hash, eqstr>::const_iterator it = reads->begin();
//					 it != reads->end(); ++it) {
//
//		const char* key = it->first;
//		printf("%s\n", key, key);
//		read_vec* read_v = it->second;
//	}
}
*/
