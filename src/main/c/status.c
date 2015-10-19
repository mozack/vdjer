#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

void print_file(char* filename, char* file_info, char* desc) {
	FILE* fp = fopen(filename, "r");
	char buf[5124];

	while (fgets(buf, 5124, fp) != NULL) {
		fprintf(stderr, "%s\t%s\t%s", file_info, desc, buf);
	}

	fflush(stdout);

	fclose(fp);
}

void print_status(char* desc) {

	print_file("/proc/self/status", "PROC_STATUS", desc);
	print_file("/proc/buddyinfo", "BUDDY_INFO", desc);
	char prefix[256];
	sprintf(prefix, "STAT:%ld", time(NULL));
	print_file("/proc/self/stat", "STAT", desc);
}
