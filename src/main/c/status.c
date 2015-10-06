#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void print_status(char* desc) {
	FILE* fp = fopen("/proc/self/status", "r");
	char buf[5124];

	while (fgets(buf, 5124, fp) != NULL) {
		fprintf(stderr, "PROC_STATUS\t%s\t%s", desc, buf);
	}

	fflush(stdout);

	fclose(fp);
}
