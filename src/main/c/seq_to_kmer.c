#include <stdio.h>
#include <stdlib.h>
#include "seq_dist.h"


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
			printf("Error converting base: %c, %d\n", base, base);
			//exit(-1);
			return -1;
	}

	return val;
}


unsigned long seq_to_int(const char* seq) {
	unsigned long val = 0;

	for (int i=0; i<SEQ_LEN; i++) {
		val = val << 2;
		int b = base_val(seq[i]);
		if (b < 0) {
			printf("seq: %s\n", seq);
			return 0;
		}
		val += b;
	}

	return val;
}

/*
int main(int argc, char** argv) {
//	ACACGGCCGTATATTT
//	CACGGCCGTATATTTC
//	CGGCCGTATATTTCTG
//	GACACGGCCGTATATT

	printf("%lu\n", seq_to_int("ACTTCTGGGGCCAGGG"));
	printf("%lu\n", seq_to_int("GACTTCTGGGGCCAGG"));
	printf("%lu\n", seq_to_int("TTCTGGGGCCAGGGAA"));
}

*/
