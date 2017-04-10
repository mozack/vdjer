import sys

id = ''
cdr3 = ''
for line in sys.stdin:
    line = line.rstrip()
    
    if line.startswith('>'):
        print line
    else:
        print line[-48:]
