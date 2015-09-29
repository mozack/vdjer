#ifndef __QUICK_MAP3__
#define __QUICK_MAP3__

#include <vector>

struct read_info {
	char* id;
	char* seq;
	char* quals;
	char read_num;
	char is_rc;
};

struct map_info {
	read_info* info;
	short pos;
};

struct read_vec {
	std::vector<read_info*>* reads;
	char* seq;
};

struct mapped_pair {
	char* contig_id;
	read_info* r1;
	read_info* r2;
	short pos1;    // read 1 position within contig
	short pos2;    // read 2 position within contig
	short insert;  // insert length
};


#endif // __QUICK_MAP3__
