import sys

sequence_fasta = sys.argv[1]
isotypes_file = sys.argv[2]
inventory = sys.argv[3]
results_dir = sys.argv[4]
vquest1 = sys.argv[5]
chain_type = sys.argv[6]
read_count_file = sys.argv[7]

def get_vq_gene(text):
    gene = 'N/A'
    genes = {}
    for field in text.split():
        if field.startswith('TRB') or field.startswith('IGH') or field.startswith('IGK') or field.startswith('IGL'):
            fields = field.split('*')
            field0 = fields[0]
            
            # Remove D from gene name
            field0 = field0.replace('D','')
            
            # Trim second - from gene name.  i.e. IGHV3-30-5 becomes IGHV3-30 
            gene_parts = field0.split('-')
            if len(gene_parts) >= 2:
                field0 = gene_parts[0] + '-' + gene_parts[1]
            
            if gene == 'N/A' or gene == field0:
                gene = field0
            elif not field0 in genes:
#                gene = 'N/A'
                gene += ',' + field0
            
            genes[field0] = 1
    
    return gene



# print header
print '\t'.join(['sample', 'sequence', 'cdr3', 'expected_counts', 'seq_id', 'isotype', 'vregion_identity', 'aa_cdr3', 'vgene', 'jgene', 'total_count'])

#
# load sequences and seq_id's
seq_ids = {}

sf = open(sequence_fasta, 'r')
id = ''
for line in sf:
    line = line.rstrip()
    if line.startswith('>'):
        id = line[1:]
    else:
        seq_ids[line] = id
        
sf.close()

#
# load sample read counts
read_counts = {}
rcf = open(read_count_file, 'r')
for line in rcf:
    line = line.rstrip()
    fields = line.split()
    if len(fields) > 1:
        sample = fields[0]
        count = int(fields[1])
        read_counts[sample] = str(count)

#
# load isotypes
isotypes = {}

if isotypes_file != "no":
    iso = open(isotypes_file, 'r')
    for line in iso:
        line = line.rstrip()
        fields = line.split()
        seq_id = fields[0]
        isotype = fields[1]
        isotypes[seq_id] = isotype
    iso.close()

#
# load vregion identity (i.e. inverse mutation load)
# and AA cdr3 from vquest
vid = {}
vq_cdr3 = {}
vq_v = {}
vq_j = {}
vidf = open(vquest1, 'r')
for line in vidf:
    if not line.startswith('Sequence'):
        line = line.rstrip()
        fields = line.split('\t')
        seq_id = fields[1]
        res = fields[2]
        vcdr3 = 'N/A'
        if res == 'No results':
            vreg_identity = 'N/A'
            vcdr3 = 'N/A'
        else:
            vreg_identity = fields[5]
            if len(fields[20].split()) > 0:
                vcdr3 = fields[20].split()[0]
        vid[seq_id] = vreg_identity
        vq_cdr3[seq_id] = vcdr3
        if (len(fields) >= 10):
            vgene = get_vq_gene(fields[3])
            jgene = get_vq_gene(fields[9])
            vq_v[seq_id] = vgene
            vq_j[seq_id] = jgene
        else:
            vq_v[seq_id] = 'N/A'
            vq_j[seq_id] = 'N/A'

#
# Go through inventory
inv = open(inventory, 'r')
for line in inv:
    sample = line.rstrip()
    # Load contigs
    sample_contigs = {}
    contig_file = results_dir + '/' + sample + '/' + chain_type + '/vdjer/vdj_contigs.fa'
    cid = ''
    cf = open(contig_file, 'r')
    for cline in cf:
        cline = cline.rstrip()
        if cline.startswith('>'):
            cid = cline[1:]
        else:
            sample_contigs[cid] = cline
    cf.close
    
    # Parse RSEM results
    rsem = results_dir + '/' + sample + '/' + chain_type + '/vdjer/rsem_results.isoforms.results'
    rf = open(rsem, 'r')
    for rline in rf:
        rline = rline.strip()
        
        if not 'transcript' in rline:
            fields = rline.split()
            cid = fields[0]
            ecount = fields[4]
            
            cid_fields = cid.split('_')
            cdr3 = cid_fields[2]
            
            # Output entry
            # sample, sequence, cdr3, expected_counts, seq_id, isotype
            if float(ecount) >= 1:
                sequence = sample_contigs[cid]
                                
                seq_id = seq_ids[sequence]
                isotype = 'N/A'
                if seq_id in isotypes:
                    isotype = isotypes[seq_id]

                if vid[seq_id] != 'N/A' and vq_cdr3[seq_id] != 'N/A':
                    output = '\t'.join([sample, sequence, cdr3, ecount, seq_id, isotype, vid[seq_id], vq_cdr3[seq_id], vq_v[seq_id], vq_j[seq_id], read_counts[sample]])
                
                print output
    rf.close()

