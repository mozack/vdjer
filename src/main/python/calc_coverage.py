#
# Calculates coverage for an input BAM file
# Includes coverage for read 1 versus read 2
# Reads must be grouped by mapped contig
# 
import sys

contig_len = int(sys.argv[1])

def output(contig, contig_len, coverage1, coverage2):
    if contig != '':
        for i in range(1,contig_len+1):
            print '\t'.join([contig, str(i), str(coverage1[i]), str(coverage2[i]), str(coverage1[i] + coverage2[i])])

def init_coverage(coverage, contig_len):
    for i in range(1,contig_len+1):
        coverage[i] = 0
        
def update_coverage(coverage, position, read_len):
    for i in range(position, position+read_len):
        coverage[i] += 1

curr_contig = ''

coverage1 = {}
coverage2 = {}

for line in sys.stdin:
    line = line.rstrip()
    fields = line.split()
    flags = int(fields[1])
    contig = fields[2]
    pos = int(fields[3])
    mate_pos = int(fields[7])
    bases = fields[9]
    read_len = len(bases)
    
    if contig != curr_contig:
        output(curr_contig, contig_len, coverage1, coverage2)
        curr_contig = contig

        init_coverage(coverage1, contig_len)
        init_coverage(coverage2, contig_len)
        
    if pos < mate_pos:
        update_coverage(coverage1, pos, read_len)
    else:
        update_coverage(coverage2, pos, read_len)
        
# Output last contig
output(contig, contig_len, coverage1, coverage2)