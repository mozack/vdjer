#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <sparsehash/sparse_hash_map>

#include "hash_utils.h"

using namespace std;

using google::sparse_hash_map;

//int READ_LEN = 50;
//int MIN_INSERT = 180 - 60; // 120
//int MAX_INSERT = 180 + 60; // 240

//TODO: Do not hardcode
int READ_LEN = 50;
int MIN_INSERT = 175;
int MAX_INSERT = 175;

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

struct read_info {
	char* id;
	char* seq;
	char* quals;
	char read_num;
	char is_rc;
};

struct read_vec {
	vector<read_info*>* reads;
	char* seq;
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

/*
void load_reads(char* file1, char* file2) {

	FILE* fp1 = fopen(file1, "r");
	FILE* fp2 = fopen(file2, "r");

	read_buf = (char*) calloc(READ_BLOCK, sizeof(char));
	read_buf_start = read_buf;

	qual_buf = (char*) calloc(READ_BLOCK, sizeof(char));
	qual_buf_start = qual_buf;

	char line1[MAX_LINE];
	char line2[MAX_LINE];
	char bases1[MAX_LINE];
	char bases2[MAX_LINE];
	char rc1[MAX_LINE];
	char rc2[MAX_LINE];
	char rc_quals1[MAX_LINE];
	char rc_quals2[MAX_LINE];

	char* read_id;
	char line_num = 1;

	while (fgets(line1, MAX_LINE, fp1) != NULL && fgets(line2, MAX_LINE, fp2) != NULL) {
		// Swallow newlines
		line1[strlen(line1)-1] = '\0';
		line2[strlen(line2)-1] = '\0';

		if (line_num == 1) {

			// Strip forward slash and anything following
			char* pch = strchr(line1, '/');
			if (pch != NULL) {
				*pch= '\0';
			}

			pch = strchr(line2, '/');
			if (pch != NULL) {
				*pch = '\0';
			}

			// Strip spaces and anything following
			pch = strchr(line1, ' ');
			if (pch != NULL) {
				*pch= '\0';
			}

			pch = strchr(line2, ' ');
			if (pch != NULL) {
				*pch = '\0';
			}

			if (strncmp(line1, line2, MAX_LINE) != 0) {
				fprintf(stderr, "Mismatched reads: %s -- %s\n", line1, line2);
				exit(-1);
			}

			// Copy read id to read buffer
			strncpy(read_buf, line1, MAX_LINE);
			read_id = read_buf;
			advance_read_buf();
		} else if (line_num == 2) {
			rc(line1, rc1);
			rc(line2, rc2);
			strncpy(bases1, line1, MAX_LINE);
			strncpy(bases2, line2, MAX_LINE);

//			printf("line1: %s -- line2: %s\n", line1, line2);
//			printf("rc1: %s -- rc2: %s\n", rc1, rc2);


		} else if (line_num == 4) {
			reverse(line1, rc_quals1);
			reverse(line2, rc_quals2);

			// Read 1, Forward
			add_read_info(read_id, bases1, line1, 1, 0);
			// Read 1, Reverse
			add_read_info(read_id, rc1, rc_quals1, 1, 1);
			// Read 2, Forward
			add_read_info(read_id, bases2, line2, 2, 0);
			// Read 2, Reverse
			add_read_info(read_id, rc2, rc_quals2, 2, 1);
		}

		if (line_num == 4) {
			line_num = 1;
		} else {
			line_num += 1;
		}
	}

	fprintf(stderr, "Num index entries: %d\n", reads->size());
	fclose(fp1);
	fclose(fp2);
}
*/

struct map_info {
	read_info* info;
	int pos;
};

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

	const char* format = "MAPPING:%s\t%d\t%s\t%d\t255\t%dM\t=\t%d\t%d\t%s\t%s\n";

	strncpy(seq, r1->info->seq, READ_LEN);
	strncpy(quals, r1->info->quals, READ_LEN);
	printf(format, read_id, flag1, contig_id, r1->pos, READ_LEN, r2->pos, insert, seq, quals);

	strncpy(seq, r2->info->seq, READ_LEN);
	strncpy(quals, r2->info->quals, READ_LEN);
	printf(format, read_id, flag2, contig_id, r2->pos, READ_LEN, r1->pos, insert, seq, quals);
}

void quick_map_process_contig(char* contig_id, char* contig) {

	int read1_count = 0;
	int MAX_READ_PAIRS = 1000000;
	map_info** read1 = (map_info**) calloc(MAX_READ_PAIRS, sizeof(map_info*));

	sparse_hash_map<const char*, struct map_info*, vjf_hash, vjf_eqstr> read2;

	// Load read 1 matches into vector
	// Load read 2 matches into map
	for (int i=0; i<strlen(contig)-READ_LEN; i++) {
		read_vec* read_v = (*reads)[contig+i];
		if (read_v != NULL) {
			for (vector<read_info*>::iterator it = read_v->reads->begin(); it != read_v->reads->end(); ++it) {
				read_info* r_info = *it;

				map_info* m_info = (map_info*) calloc(1, sizeof(map_info));
				m_info->info = r_info;
				m_info->pos = i + 1;

				if (r_info->read_num == 1) {
//					printf("Adding R1: %s\n", r_info->id);
					read1[read1_count++] = m_info;
				} else {
//					printf("Adding R2: %s\n", r_info->id);
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
				int insert = abs(r1->pos - r2->pos) + READ_LEN;
				if (insert >= MIN_INSERT && insert <= MAX_INSERT) {
					// We have a hit.  Output
					output_mapping(contig_id, r1, r2, insert);
				}
			}
		}
	}

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

void process_contigs(char* contig_file) {
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

			quick_map_process_contig(contig_id, contig);
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

void output_header(char* contig_file, int argc, char** argv) {
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


	char pg[10000];
	pg[0] = '\0';
	int chars_left = 9999;
	for (int i=0; i<argc; i++) {
		strncat(pg, argv[i], chars_left);
		strcat(pg, " ");
		chars_left = 10000 - strlen(pg);
	}

	printf("@PG\tID:quickmap\tPN:quickmap\tCL:%s\n", pg);

	fclose(fp);
	fprintf(stderr, "Done writing header\n");
	fflush(stderr);
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
