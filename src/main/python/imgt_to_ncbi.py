# Joins an IMGT gene list with an NCBI GFF.
# IMGT genes downloaded from here: http://www.imgt.org/genedb/doc#directlinksgroup
# NCBI GFF downloaded from here: ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/ARCHIVE/BUILD.37.3/GFF/ref_GRCh37.p5_top_level.gff3.gz

import sys

imgt = sys.argv[1]
ncbi = sys.argv[2]

#imgt = '/home/lmose/dev/vdj/genes/ig.f.txt'
#ncbi = '/home/lmose/dev/vdj/genes/ref_GRCh37.p5_top_level.gff3'

#imgt = '/home/lmose/dev/vdj/genes/28454.test'
#ncbi = '/home/lmose/dev/vdj/genes/28454.ncbi.test'

def get_chr(input):
    if not input.startswith('NC'):
        return None
    else:
        start = input.find('0')
        stop = input.find('.')
        chr = 'chr' + str(int(input[start:stop]))
        return chr

IMGT = open(imgt, 'r')

imgt_recs = {}

for line in IMGT:
    line = line.rstrip()
    if not line.startswith('Species'):
        imgt_fields = line.split(';')
        ncbi_id = imgt_fields[9]
        imgt_recs[ncbi_id] = imgt_fields
    
IMGT.close()

NCBI = open(ncbi, 'r')
for line in NCBI:
    line = line.rstrip()
    if not line.startswith('#'):
        ncbi_fields = line.split('\t')
        attributes = ncbi_fields[8]
        idx = attributes.find('Dbxref=GeneID:')
        if idx > 0 and ncbi_fields[2] == 'gene':
            idx += len('Dbxref=GeneID:')
            stop = attributes[idx:].find(',')
            stop2 = attributes[idx:].find(';')
            if stop2 > 0 and (stop2 < stop or stop < 0):
                stop = stop2
            if stop > 0:
                id = attributes[idx:stop+idx]
                
                if not id == None and not id == '' and id in imgt_recs:
                    rec = imgt_recs[id]
                    
                    chr = get_chr(ncbi_fields[0])
                    
                    if chr:
                        print '\t'.join([chr, ncbi_fields[3], ncbi_fields[4], rec[1], '.', ncbi_fields[6], rec[2], ncbi_fields[2], id])

NCBI.close()
    