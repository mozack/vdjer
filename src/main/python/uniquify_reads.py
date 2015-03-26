import sys

count = 1

for line in sys.stdin:
    line = line.rstrip()
    if line.startswith('@'):
        print line
    else:
        fields = line.split('\t')
        fields[0] = fields[0] + '_' + str(count)
        count += 1
        print '\t'.join(fields)
    