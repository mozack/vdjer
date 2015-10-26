#ifndef __HASH_UTILS__
#define __HASH_UTILS__

extern int kmer_size;
extern int CONTIG_SIZE;
extern int VREGION_KMER_SIZE;

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


//
// Equals operator and hash on string using global "CONTIG_SIZE" variable.
//
struct contig_eqstr
{
  bool operator()(const char* s1, const char* s2) const
  {
    return (s1 == s2) || (s1 && s2 && strncmp(s1, s2, CONTIG_SIZE) == 0);
  }
};

struct contig_hash
{
	uint64_t operator()(const char* kmer) const
	{
		return MurmurHash64A(kmer, CONTIG_SIZE, 97);
		//return chunk;
	}
};


//
// Equals operator and hash on string using global "VREGION_KMER_SIZE" variable.
//
struct vregion_eqstr
{
  bool operator()(const char* s1, const char* s2) const
  {
    return (s1 == s2) || (s1 && s2 && strncmp(s1, s2, VREGION_KMER_SIZE) == 0);
  }
};

struct vregion_hash
{
	uint64_t operator()(const char* kmer) const
	{
		return MurmurHash64A(kmer, VREGION_KMER_SIZE, 97);
		//return chunk;
	}
};

#endif // __HASH_UTILS__

