#
# Extracts anchor kmers from a fasta file supplied as stdin
# usage cat my.fa | python anchor_kmers.py <side> <kmer> <offset>
# side = 5 or 3
# kmer = kmer len
# offet = offset from 5' or 3' end

import sys

side = int(sys.argv[1])
kmer_len = int(sys.argv[2])
offset = int(sys.argv[3])

start = 0
stop = 0

if side == 3:
    start = -(kmer_len+offset)
    stop  = -offset
elif side == 5:
    start = offset
    stop  = kmer_len+offset
else:
    print 'bad side: ' + side

for line in sys.stdin:
    line = line.rstrip()
    
    if not line.startswith('>'):
        kmer = line[start:stop]
        print kmer