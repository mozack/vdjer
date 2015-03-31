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
