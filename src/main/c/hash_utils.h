#ifndef __HASH_UTILS__
#define __HASH_UTILS__

extern int kmer_size;

uint64_t MurmurHash64A ( const void * key, int len, uint64_t seed );

//
// Equals operator and hash on string using global "kmer_size" variable.
//
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

//
// Equals operator and hash on string using standard string comparison
//
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

#endif // __HASH_UTILS__
