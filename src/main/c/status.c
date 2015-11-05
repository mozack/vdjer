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

time_t status_start_time = time(NULL);
time_t status_prev_time = time(NULL);

void print_status(char* desc) {

	time_t curr_time = time(NULL);
	fprintf(stderr, "ELAPSED_SECS\t%s\t%ld\t%ld\n", desc, curr_time - status_start_time, curr_time-status_prev_time);
	status_prev_time = curr_time;
	print_file("/proc/self/status", "PROC_STATUS", desc);
	print_file("/proc/buddyinfo", "BUDDY_INFO", desc);
	char prefix[256];
	sprintf(prefix, "STAT:%ld", time(NULL));
	print_file("/proc/self/stat", prefix, desc);
}
