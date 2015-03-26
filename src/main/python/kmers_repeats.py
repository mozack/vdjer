import sys

MIN_BASE_QUAL = 13
KMER = 15

ctr = 0
bases = ''
quals = ''
id = ''
sep = ''


for line in sys.stdin:
    kmers = {}
    line = line.rstrip()
    ctr += 1
    if ctr == 1:
        id = line
    elif ctr == 2:
        bases = line
    elif ctr == 3:
        sep = line
    elif ctr == 4:
        ctr = 0
        quals = line
        
        for i in range(len(bases)-KMER):
            skip = False
            kmer = bases[i:i+KMER]
            
            kmer_quals = quals[i:i+KMER]
            for q in kmer_quals:
                if (ord(q)-ord('!')) < MIN_BASE_QUAL:
                    skip = True

            if not skip:            
                if kmer in kmers:
                    kmers[kmer] = kmers[kmer]+1
                    print '\t'.join(['REPEAT',kmer,str(kmers[kmer]),bases])
                else:
                    kmers[kmer] = 1
 