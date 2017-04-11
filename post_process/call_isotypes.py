import sys

bed=sys.argv[1]

chromosome = 'chr14'

lookup = {}
results = {}

b=open(bed, 'r')
for line in b:
    line = line.rstrip()
    fields = line.split();
    start = int(fields[1])
    stop = int(fields[2])
    isotype = fields[3]
    for i in range(start,stop):
        lookup[i] = isotype
        
for line in sys.stdin:
    line = line.rstrip()
    fields = line.split()
    id = fields[0]
    pos = int(fields[3])
    
    if pos in lookup:
        if id not in results:
            # Add new results
            results[id] = [ lookup[pos] ]
        else:
            results[id].append(lookup[pos])
            
for id in results:
    call = ','.join(results[id])
    type = {}
    for text in results[id]:
        subtype = text[0:4]
        type[subtype] = 1
                
    print id + '\t' + ','.join(type.keys()) + '\t' + call

