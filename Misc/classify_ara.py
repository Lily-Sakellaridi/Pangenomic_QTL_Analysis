"""
Script to classify Arabidopsis homology groups into core, core & single-copy orthologs,
accessory, unique, or non-variant.
Author: Lily Sakellaridi.
"""

from sys import argv

def make_group_dict(groups_file):
    """
    Make a dictionary where key: groups, value: list of gene counts per genome per group.    
    """
    group_dict = {}
    with open(groups_file) as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip().split(": ")
            group_id = line[0]
            genomes = [0] * 19 # value
            #print(genomes)
            #print(group_id)
            genes = line[1].split(" ")
            for gene in genes:
                gene = gene.split("#")
                #print(gene)
                #gene_id = gene[0]
                genome = int(gene[1])
                #print(genome)
                index = genome - 1 # compensate for Python's zero-based indexing
                genomes[index] += 1
            #print(genomes)
            group_dict[group_id] = genomes
            
    return group_dict
    
def classify(group_dict):
    """
    Classify groups into core, core & sco, accessory and unique.
    """
    new_group_dict = {}
    for group, genomes in group_dict.items():
        type = "Undefined type" # Initialize type; if you actually get 'undefined type', something is wrong.
        variant_type = "Undefined variation type"
        if 0 in genomes: # This is either a presence/absence or a complex event.
            if genomes.count(0) == 18: # a unique group cannot have copy number variation.
                type = "unique"
                variation_type = "pav"
            else:
                type = "accessory"
                values = set(genomes)
                if len(values) > 2:  # this means there is both presence-absence and copy number variation
                    variation_type = "complex"
                else:
                    assert len(values) == 2 # this means there is only presence-absence variation.
                    variation_type = "pav"
        else:
            values = set(genomes)
            if len(values) == 1:
                variation_type = "non-variant"
                if 1 in values:
                    type = "core & sco"
                else:
                    type = "core"
            else:
                type = "core"
                variation_type = "cnv"
                
        new_group_dict[group] = [genomes, type, variation_type]
        
    return new_group_dict

if __name__=="__main__":
    groups = argv[1]
    gd = make_group_dict(groups)
    #print(gd)
        
    ngd = classify(gd)
    for k, v in ngd.items():
        print(k, v)