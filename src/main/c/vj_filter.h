#ifndef __VJ_FILTER__
#define __VJ_FILTER__

#define BIG_CONSTANT(x) (x##LLU)

uint64_t MurmurHash64A ( const void * key, int len, uint64_t seed );

struct vjf_eqstr
{
  bool operator()(const char* s1, const char* s2) const
  {
    return (s1 == s2) || (s1 && s2 && strcmp(s1, s2) == 0);
  }
};

struct vjf_hash
{
	uint64_t operator()(const char* seq) const
	{
		uint64_t h = MurmurHash64A(seq, strlen(seq), 97);
		return h;
	}
};

// Init V and J anchor indices as well as search params
void vjf_init(char* v_file, char* j_file, int max_dist, int min_win, int max_win,
		char j_conserved, int window_span, int j_extension);

// Search for candidate VDJ windows
void vjf_search(char* contig, google::sparse_hash_set<const char*, vjf_hash, vjf_eqstr>& windows);

#endif
