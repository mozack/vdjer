import sys
import subprocess

KMER=25

#VDJG='/datastore/nextgenout4/seqware-analysis/lmose/vdj/mouse/abra/vdjg'
VDJG='/datastore/nextgenout4/seqware-analysis/lmose/vdj/mouse/samples/vdjg.bash'

rc_bases = {}
rc_bases['A'] = 'T'
rc_bases['T'] = 'A'
rc_bases['C'] = 'G'
rc_bases['G'] = 'C'

aa_bases = {}

def add(aa_bases, aa, bases_str):
    bases_list = bases_str.split(',')
    for bases in bases_list:
        aa_bases[bases] = aa

add(aa_bases,'A','GCA,GCC,GCG,GCT')
add(aa_bases,'C','TGC,TGT')
add(aa_bases,'D','GAC,GAT')
add(aa_bases,'E','GAA,GAG')
add(aa_bases,'F','TTC,TTT')
add(aa_bases,'G','GGA,GGC,GGG,GGT')
add(aa_bases,'H','CAC,CAT')
add(aa_bases,'I','ATA,ATC,ATT')
add(aa_bases,'K','AAA,AAG')
add(aa_bases,'L','TTA,TTG,CTA,CTC,CTG,CTT')
add(aa_bases,'M','ATG')
add(aa_bases,'N','AAC,AAT')
add(aa_bases,'P','CCA,CCC,CCG,CCT')
add(aa_bases,'Q','CAA,CAG')
add(aa_bases,'R','CGA,CGC,CGG,CGT,AGA,AGG')
add(aa_bases,'S','TCA,TCC,TCG,TCT,AGC,AGT')
add(aa_bases,'T','ACA,ACC,ACG,ACT')
add(aa_bases,'V','GTA,GTC,GTG,GTT')
add(aa_bases,'W','TGG')
add(aa_bases,'Y','TAC,TAT')
add(aa_bases,'.','TAA,TAG,TGA')


# Convert input bases to Amino Acids
def convert_to_aa(bases):
    aminos = ''
    for i in range(0,len(bases)-2, 3):
        aminos += aa_bases[bases[i:i+3]]
        
    return aminos
    

# Reverse Complement
def rc(bases):
    bases = bases[::-1]
    new_bases = ''
    
    for base in bases:
        new_bases += rc_bases[base]
    
    return new_bases

# There should be 1 line in format <score>    <window(aminos)>    <window(nucleotides)>    <window(nucleotides)>    ...
for line in sys.stdin:
    line = line.rstrip()
    fields = line.split('\t')
        
    # For each window
    for i in range(2,len(fields)):
        window = fields[i]

        # Get kmer for beginning and end of window
        start_kmer = window[0:KMER]
        end_kmer = window[-KMER:]
        start_kmer = rc(start_kmer)

        start = ''
        stop = ''
        
        # Invoke assembler script passing in start kmer and capture stdout
        # Extends 5' side of window
        vdj = subprocess.check_output([VDJG, start_kmer])
        vdj = vdj.rstrip()
        
        # Get 4th column from vdj
        vdj_fields = vdj.split('\t')
        if (len(vdj_fields) >= 4):
            # Strip off leading kmer and reverse complement
            start = rc(vdj_fields[3][25:])
            
        # Invoke assembler script passing in end kmer and capture stdout
        # Extends 3' side of window
        vdj = subprocess.check_output([VDJG, end_kmer])
        vdj = vdj.rstrip()
        
        vdj_fields = vdj.split('\t')
        if (len(vdj_fields) >= 4):
            # Strip off leading kmer
            stop = vdj_fields[3][25:]

        if start != '' and stop != '':
            
            # Output AA and bases
            rna = start + window + stop
            amino_acids = convert_to_aa(rna)
            print amino_acids + '\t' + rna
            

