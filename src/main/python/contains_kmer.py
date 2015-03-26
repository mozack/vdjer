import sys

seq = sys.argv[1]

for line in sys.stdin:
    line = line.rstrip()
    
    fields = line.split('\t')
    
    if seq.find[fields[0]] > -1:
        print line + '\t1'
    else:
        print line + '\t0'