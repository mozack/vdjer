import sys

# Assuming first line is the top CDR
is_best = True
count = 1
for line in sys.stdin:
    line = line.rstrip()
    fields = line.split('\t')
    
    name = 'contig_';
    if is_best:
        name = 'contig_best_'
        is_best = False
    
    for i in range(2,len(fields)): 
        print '>' + name + str(count)
        print fields[i]
        count += 1 
        