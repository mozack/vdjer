# Filters IMGT gene entries that do not contribute to functional genes
# IMPORTANT: Double check function column values if applying to new datasets
# F - appears to always come first in the current dataset
import sys

for line in sys.stdin:
    line = line.rstrip()
    
    if line.startswith('Species'):
        print line
    else:
        fields = line.split(';')
        func = fields[2]
        if func.startswith('F'):
            print line