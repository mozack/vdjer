import sys

tcga_id = sys.argv[1]

total_reads = -1
cdr3_reads = 0
best_cdr3_reads = 0
window_count = 0

def get_val(file):
    f = open(file, 'r')
    line = f.readline().rstrip()
    val = float(line)
    f.close()
    
    return val

f = open('../cdr3.best.txt')
line = f.readline().rstrip()
fields = line.split('\t')
cdr3 = ''
if len(fields) >= 2:
    cdr3 = fields[1]

total_reads = get_val('total.readcount.txt')

if cdr3 != '':
    cdr3_reads = get_val('cdr3.readcount.txt')
    best_cdr3_reads = get_val('best.readcount.txt')
    window_count = get_val('window_count.txt')

cdr3_frac = cdr3_reads / total_reads
best_frac = 0
if cdr3_reads > 0:
    best_frac = best_cdr3_reads / cdr3_reads
    
best_frac2 = 0
if total_reads > 0:
    best_frac2 = best_cdr3_reads / total_reads

values = '\t'.join([tcga_id,str(window_count),str(total_reads), str(cdr3_reads),str(best_cdr3_reads),str(cdr3_frac),str(best_frac),str(best_frac2),cdr3])

print values