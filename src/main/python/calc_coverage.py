#
# Calculates coverage for an input BAM file
# Includes coverage for read 1 versus read 2
# Reads must be grouped by mapped contig
# 
import sys

contig_len = int(sys.argv[1])

def output(contig, contig_len, coverage1, coverage2, insert1, insert2):
    if contig != '':
        for i in range(1,contig_len+1):
            i1 = 0
            if coverage1[i] > 0:
                i1 = insert1[i] / coverage1[i]
            i2 = 0
            if coverage2[i] > 0:
                i2 = insert2[i] / coverage2[i]
            print '\t'.join([contig, str(i), str(coverage1[i]), str(coverage2[i]), str(coverage1[i] + coverage2[i]), str(i1), str(i2)])

def init_coverage(coverage, contig_len):
    for i in range(1,contig_len+1):
        coverage[i] = 0
        
def update_coverage(coverage, position, read_len, insert_list, insert_len):
    for i in range(position, position+read_len):
        coverage[i] += 1
        insert_list[i] += insert_len

curr_contig = ''

coverage1 = {}
coverage2 = {}

insert1 = {}
insert2 = {}

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
        output(curr_contig, contig_len, coverage1, coverage2, insert1, insert2)
        curr_contig = contig

        init_coverage(coverage1, contig_len)
        init_coverage(coverage2, contig_len)
        init_coverage(insert1, contig_len)
        init_coverage(insert2, contig_len)
        
    if pos < mate_pos:
        update_coverage(coverage1, pos, read_len, insert1, insert)
    else:
        update_coverage(coverage2, pos, read_len, insert2, insert)
        
# Output last contig
output(contig, contig_len, coverage1, coverage2, insert1, insert2)