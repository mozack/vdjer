# Remove first N bases / quality scores from fastq file 
import sys

n1 = int(sys.argv[1])
n2 = int(sys.argv[2])
file1 = sys.argv[3]
file2 = sys.argv[4]
out1 = sys.argv[5]
out2 = sys.argv[6]

f1 = open(file1, 'r')
f2 = open(file2, 'r')

o1 = open(out1, 'w')
o2 = open(out2, 'w')

count = 1

for line1 in f1:
    line1 = line1.rstrip()
    line2 = f2.readline().rstrip()
    
    if count == 2 or count == 4:
        line1 = line1[n1:]
        line2 = line2[n2:]
        if len(line2) > len(line1):
            line2 = line2[:len(line1)]

    o1.write(line1)
    o1.write('\n')
    o2.write(line2)
    o2.write('\n')
    
    if count == 4:
        count = 1
    else:
        count += 1
        
f1.close()
f2.close()
o1.close()
o2.close()