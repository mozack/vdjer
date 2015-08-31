#
# Calculates coverage for an input BAM file
# recording read start positions spanning each
# position and mate read start positions
# Reads must be grouped by mapped contig
#
# Output columns:  contig_id, position, read start positions, mate start positions

import sys

contig_len = int(sys.argv[1])

def output(contig, contig_len, read1, read1_mate, read2, read2_mate):
    if contig != '':
        for i in range(1,contig_len+1):
            r1 = ','.join(read1[i])
            m1 = ','.join(read1_mate[i])
            r2 = ','.join(read2[i])
            m2 = ','.join(read2_mate[i])
            if r1 == '':
                r1 = '.'
            if m1 == '':
                m1 = '.'
            if r2 == '':
                r2 = '.'
            if m2 == '':
                m2 = '.'
            
#            print contig + '\t' + str(i) + '\t' + r1 + '\t' + m1 + '\t' + r2 + '\t' + m2
            
            print '\t'.join([contig, str(i), r1, m1, r2, m2])

def init_coverage(coverage, contig_len):
    for i in range(1,contig_len+1):
        coverage[i] = []
        
def update_coverage(reads, read_mates, read_start, mate_start, read_len):
    for i in range(read_start, read_start+read_len):
        reads[i].append(str(read_start))
        read_mates[i].append(str(mate_start))

curr_contig = ''

read1 = {}
read1_mate = {}
read2 = {}
read2_mate = {}

for line in sys.stdin:
    line = line.rstrip()
    fields = line.split()
    flags = int(fields[1])
    contig = fields[2]
    pos = int(fields[3])
    mate_pos = int(fields[7])
    insert = abs(int(fields[8]))
    bases = fields[9]
    read_len = len(bases)
    
    if contig != curr_contig:
        output(curr_contig, contig_len, read1, read1_mate, read2, read2_mate)
        curr_contig = contig

        init_coverage(read1, contig_len)
        init_coverage(read1_mate, contig_len)
        init_coverage(read2, contig_len)
        init_coverage(read2_mate, contig_len)
        
    if pos < mate_pos:
        update_coverage(read1, read1_mate, pos, mate_pos, read_len)
    else:
        update_coverage(read2, read2_mate, pos, mate_pos, read_len)
        
# Output last contig
output(curr_contig, contig_len, read1, read1_mate, read2, read2_mate)