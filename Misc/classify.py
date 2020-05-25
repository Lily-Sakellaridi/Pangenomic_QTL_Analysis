"""
Script to: 1) make the correct classified_hgroups table on gene level
from the pantools_homology_groups file, 2) add this information
to the QTL tables.

Author: Lily Sakellaridi
"""

from sys import argv

def add_info_to_table(qtl_table, group_dict):
    """
    Add group info to the QTL table.
    """
    prefix_name = qtl_table[:-4]
    with open(qtl_table) as qf:
        lines = qf.readlines()
        with open(prefix_name + "_with_group_info", "w") as out:
            out.write(lines[0].strip() + ", type, members_g1, members_g2, members_g3\n")
            for line in lines[1:]:
                line = line.strip().split(",")
                #qtl_id = line[0]
                group_id = line[0].strip()
                info = group_dict[group_id]
                group_type = info[0]
                genome1 = str(info[1])
                genome2 = str(info[2])
                genome3 = str(info[3])
                line_end = [group_type, genome1, genome2, genome3]
                line = line + line_end
                #print(line)
                #print(group_id, group_type, genome1, genome2, genome3)
                #line = line + info
                out.write(",".join(line) + "\n")

def write_output(group_dict):
    with open("classified_groups_gene_level", "w") as out:
        for k, v in group_dict.items():
            group_id = k
            group_type = v[0]
            group_genome1 = str(v[1])
            group_genome2 = str(v[2])
            group_genome3 = str(v[3])
            row = [group_id, group_type, group_genome1, group_genome2, group_genome3]
            row = ",".join(row)
            out.write(row + "\n")
            

def classify(infile):
    """
    Return a dictionary where keys: group identifiers, values: information lists.
    
    Each list contains four items:
    1. group type (str): 'core & not sco', 'core & sco', 'accessory', 'unique'
    2. number of genes on genome 1 (int)
    3. number of genes on genome 2 (int)
    4. number of genes on genome 3 (int)
    """
    with open(infile) as f:
        lines = f.readlines()
        group_dict = {}
        for line in lines:
            genome_1 = set()
            genome_2 = set()
            genome_3 = set()
            line = line.strip().split()
            group_id = line[0][:-1]
            transcripts = line[1:]
            for transcript in transcripts:
                gene_id = transcript[:-4]
                if transcript.endswith("1"):
                    genome_1.add(gene_id)
                elif transcript.endswith("2"):
                    genome_2.add(gene_id)
                elif transcript.endswith("3"):
                    genome_3.add(gene_id)
                else:
                    print("Nonsensical transcript genome number.")
                    break
                    
            members_g1 = len(genome_1)
            members_g2 = len(genome_2)
            members_g3 = len(genome_3)
            
            group_type = infer_type(members_g1, members_g2, members_g3)
            group_dict[group_id] = (group_type, members_g1, members_g2, members_g3)
            #yield(group_id, group_type, members_g1, members_g2, members_g3)
            
    return group_dict
            
def infer_type(num_1, num_2, num_3):
    sco = ((num_1 == 1) and (num_2 == 1) and (num_3 == 1))
    accessory = (((num_1 == 0) and (num_2 != 0) and (num_3 != 0))
                 or ((num_1 != 0) and (num_2 != 0) and (num_3 == 0))
                 or ((num_1 != 0) and (num_2 == 0) and (num_3 != 0)))
    unique = (((num_1 == 0) and (num_2 == 0) and (num_3 != 0))
              or ((num_1 == 0) and (num_2 != 0) and (num_3 == 0))
              or ((num_1 != 0) and (num_2 == 0) and (num_3 == 0)))
                 
    if sco:
        return "core & sco"
    elif accessory:
        return "accessory"
    elif unique:
        return "unique"
    else:
        return "core & not sco"

if __name__=="__main__":
    hgroups = argv[1]
    qtl = argv[2]
    
    gd = classify(hgroups)
    #print(gd.keys())
    #write_output(gd)
    add_info_to_table(qtl, gd)
