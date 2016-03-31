import sys

next_cluster_id = 1
clusters = {}

for line in sys.stdin:
    line = line.rstrip()
    
    if line.startswith('sample'):
        print line + '\t' + 'cluster'
    else:
        fields = line.split('\t')
        sample = fields[0]
        isotype = fields[5]
        cdr3 = fields[7]
        vgene = fields[8]
        jgene = fields[9]
        
        cluster_key = '_'.join([sample, isotype, cdr3, vgene, jgene])
        
        if cluster_key not in clusters:
            clusters[cluster_key] = next_cluster_id
            next_cluster_id += 1
        
        print line + '\tcls_' + str(clusters[cluster_key])