import sys

KMER = 150

ctr = 0
bases = ''
id = ''

kmers = {}

included_kmers = {}
restrict_kmers = False

if (len(sys.argv) > 1):

    restrict_kmers = True
    f = open(sys.argv[1],'r')
#    f = open('/home/lmose/dev/vdj/kmer_test/freq.txt','r')
    for line in f:
        line = line.rstrip()
        fields = line.split('\t')
        if len(fields) > 0:
            kmer = fields[0]
            included_kmers[kmer] = 1
    

#f2 = open('/home/lmose/dev/vdj/kmer_test/test.fastq', 'r')
#for line in f2:

for line in sys.stdin:
    line = line.rstrip()
    ctr += 1
    if ctr == 1:
        id = line
    elif ctr == 2:
        bases = line
        ctr = 0
        
        for i in range(len(bases)-KMER):
            skip = False
            kmer = bases[i:i+KMER]
            
            if restrict_kmers:
                if not kmer in included_kmers:
                    skip = True
                                
            if not skip:
                if kmer in kmers:
                    kmers[kmer] = kmers[kmer]+1
                else:
                    kmers[kmer] = 1

if len(included_kmers) > 0:
    for kmer in included_kmers:
        if kmer in kmers:
            print kmer + '\t' + str(kmers[kmer])
        else:
            print kmer + '\t0'
else:
    for kmer in kmers:
        if kmers[kmer] >= 1000:
            print kmer + '\t' + str(kmers[kmer])
     