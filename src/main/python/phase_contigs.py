import sys

READ_LEN = 50
INSERT_25 = 163
INSERT_75 = 197
READ_SPAN = 40
READ_LEFT = READ_LEN-READ_SPAN

EVAL_START = 25
EVAL_STOP = 390

FLOOR = 2

curr_contig = ''
is_filtered = False

for line in sys.stdin:
    line = line.rstrip()
    fields = line.split('\t')
    contig = fields[0]
    pos = int(fields[1])
    
    if contig != curr_contig:
        if not is_filtered and curr_contig != '':
            print '\t'.join([curr_contig, 'PASSED'])
        curr_contig = contig
        is_filtered = False
    
    if is_filtered or pos < EVAL_START or pos > EVAL_STOP:
        continue
    
    if len(fields) < 5:
        print '\t'.join([contig, 'FILTER', 'READ_COV', str(pos)])
        is_filtered = True
        continue
    
    # Phase reads 
    reads = []
    if fields[2] != '.':
        reads.extend(map(int, fields[2].split(',')))
    
#    print 'reads1: ' + '\t' + str(reads)
    
    if fields[4] != '.':
        reads.extend(map(int, fields[4].split(',')))
    
#    print 'reads1: ' + '\t' + str(reads)
    
    read_floor = 0
    read_target = READ_LEN - READ_SPAN

    # i.e. 25 - (390-40)
    if pos < EVAL_STOP - READ_SPAN:
        for read in reads:
             if pos - read < READ_LEFT:
                 read_floor += 1
                 if read_floor >= FLOOR:
                     break
    
        # Filter if read floor not reached
        if read_floor < FLOOR:
            print '\t'.join([contig, 'FILTER', 'READ_FLOOR', str(pos)])
            is_filtered = True
            continue
    
    # Phase mates
    
    
    
#    mates = map(int, fields[3].split(','))
    
if not is_filtered and curr_contig != '':
    print '\t'.join([contig, 'PASSED'])
