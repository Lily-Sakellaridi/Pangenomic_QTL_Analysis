"""
===========
Visualize QTL and their homologous regions, including
links between individual genes.
===========

Chromosomes and scaffolds are plotted as a broken bar collection.
Genes are plotted as points on the bars.
Genes in the same group are linked with lines.
Prioritized genes from strigolactone use case are plotted as cyan lines.
Links to genes on other chromosomes are not shown.

Author: Lily Sakellaridi
"""

from sys import argv
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

def parse_infile(infile, start, y, direction, scaffold_end, is_second, real_start):
    with open(infile) as f:
        lines = f.readlines()
        for line in lines:
            gene = line.strip().split(",")
            #print(start)
            x = int(gene[3]) - start
            
            if direction == "reverse":
                if is_second == "yes":
                    start = real_start + 50000
                    x = int(gene[3]) - start
                print("start:" + str(start))
                print("scaffold end: " + str(scaffold_end))
                print("x coordinate before: " + str(x))
                x = scaffold_end - x
                print("x coordinate after: " + str(x))
            x += 10000 #push genes a bit to the right
            gene = [gene[0], int(gene[2]), x, y, gene[5], gene[1]]
            yield gene
            
def get_prioritized_genes(pri_file):
    prio_genes = []
    with open(pri_file) as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip().split(" ")
            prio_genes.append(line[3])
            
    return prio_genes
            
def plot_genes(genes, group_dict):
    for gene in genes:
        #group = gene[0]
        #sc = gene[1]
        x = gene[2]
        y = gene[3]
        locus = gene[5]
        if locus in prioritized_genes:
            #plt.plot(x, y - 3, "cv")
            plt.plot(x, y, "c|", markersize = 38)
        #is_transposon = gene[4]
        #marker_color = set_marker_shape(gene)
        marker = set_marker(gene, group_dict)
        #marker = marker_color + marker_shape
        plt.plot(x, y, marker, alpha = 0.8)
        #if 'TRUE' in is_transposon:
        #    plt.plot(x, y, 'r*')
        #else:
        #    plt.plot(x, y, 'g*')
            
def make_group_dict(groups_file):
    group_dict = {}
    with open(groups_file) as g:
        lines = g.readlines()
        for line in lines:
            line = line.strip().split(",")
            group = line[0]
            type = line[1]
            group_dict[group] = type
            
    return group_dict
    
def plot_links(group1, group2):
    for point1 in group1:
        for point2 in group2:
            if point1[0] == point2[0]:  # point[0] is the homology group
                #print(point1[0], point2[0])
                
                plt.plot([point1[2], point2[2]], [point1[3], point2[3]], 'k-', alpha=0.2)
            
def set_marker(gene, group_dict):
    """
    Set marker color and shape according to gene type.
    
    Red square: core and single copy ortholog.
    Blue circle: core (and not single copy ortholog).
    Green thin diamond: accessory.
    Magenta star: unique.
    """
    group = gene[0]
    type = group_dict[group]
    #print(type)
    
    if "core & sco" in type:
        return "r+"
    if "core & not sco" in type:
        return "bx"
    elif "accessory" in type:
        return "g."
    elif "unique" in type:
        return "m*"
    else:
        return "k" # If black markers appear, there is something wrong, since there are no other types.
        
        
if __name__=="__main__":
    
    # Take in infiles.
    nipponbare = "Nipponbare"
    dj123_f1 = "dj123_f1.csv" 
    ir64_f1 = "ir64_f1" 
    dj123_f2 = "dj123_f2.csv" 
    ir64_f2 = "ir64_f2" 
    classified_groups = "classified_groups_gene_level"
    prioritized_genes = get_prioritized_genes("prioritized_genes.txt")

    fig, ax = plt.subplots()
    
    nipponbare_bottom = 30 # Nipponbare in the middle
    dj123_bottom = 50
    ir64_bottom = 10
    
    bar_width = 6
    
    # (Real) starts of regions corresponding to first bars.
    # Used in the parse_infile function to adjust the x coordinate of the gene.
    # (Starts on the plot are 0.)
    nipponbare_start = 29000000
    dj123_f1_start = 2676
    dj123_f2_start = 248610
    ir64_f1_start = 163080
    ir64_f2_start = 12975
    
    
    # Sizes of QTL and homologous regions.
    qtl_size = 620000 # on Nipponbare
    dj123_f1_size = 238769 - dj123_f1_start + 20000
    dj123_f2_size =  525316 - dj123_f2_start + 20000
    ir64_f1_size =  386801 - ir64_f1_start + 20000
    ir64_f2_size =  295213 - ir64_f2_start + 20000

    # Distance between fragments. Arbitrary number, must be consistent.
    dist = 50000
    
    # Starts (on the plot) of regions of subsequent bars.
    # Used in the parse_infile function to set the y coordinate of the gene.
    # (Starts of first bars are 0.)
    dj123_f2_plot_start = dj123_f1_size + dist # sc 2
    ir64_f2_plot_start = ir64_f1_size + dist # sc 142
    
    # Parse group classification file.
    group_dict = make_group_dict(classified_groups)
    
    # Parse gene files.
    
    nipponbare_genes = list(parse_infile(nipponbare, nipponbare_start, nipponbare_bottom + 3, "forward", 1, "no", 1))
    dj123_f1_genes = list(parse_infile(dj123_f1, dj123_f1_start,dj123_bottom + 3, "reverse", 238769, "no", 1))
    dj123_f2_genes = list(parse_infile(dj123_f2, dj123_f2_start - dj123_f2_plot_start,dj123_bottom + 3, "reverse", 525316, "yes", 248610))
    ir64_f1_genes = list(parse_infile(ir64_f1, ir64_f1_start, ir64_bottom + 3, "forward", 1, "no", 1))
    ir64_f2_genes = list(parse_infile(ir64_f2, ir64_f2_start - ir64_f2_plot_start, ir64_bottom + 3, "forward", 1, "no", 1))
    
    
    dj123_genes = dj123_f1_genes + dj123_f2_genes 
    ir64_genes = ir64_f1_genes + ir64_f2_genes

    # Plot bars.
    nipponbare_bar = ax.broken_barh([(0, qtl_size)], (nipponbare_bottom, bar_width), facecolors='beige')
    dj123_f1_bar = ax.broken_barh([(0, dj123_f1_size), (dj123_f2_plot_start, dj123_f2_size)], (dj123_bottom, bar_width), facecolors = 'beige')
    ir64_f1_bar = ax.broken_barh([(0, ir64_f1_size), (ir64_f2_plot_start, ir64_f2_size)], (ir64_bottom, bar_width), facecolors = 'beige')
    
    ax.set_ylim(5, 80)
    ax.set_xlim(1, 625000)
    ax.set_yticks([ir64_bottom + 3, dj123_bottom + 3, nipponbare_bottom + 3])
    ax.set_yticklabels(['IR64', 'DJ123', 'Nipponbare'])
    
    # Plot genes.
    plot_genes(nipponbare_genes, group_dict)
    plot_genes(dj123_f1_genes, group_dict) # sc end = 238769
    plot_genes(ir64_f1_genes, group_dict)
    plot_genes(dj123_f2_genes, group_dict) # sc end = 525316
    plot_genes(ir64_f2_genes, group_dict)
    
    # Plot links between genes.
    plot_links(nipponbare_genes, dj123_genes)
    plot_links(nipponbare_genes, ir64_genes)
    #plot_links(nipponbare_genes, nipponbare_genes)
    
    # Add label that explains gene types.
    blue_x = mlines.Line2D([], [], color='blue', marker='x', linestyle='None',
                          markersize=10, label='Core')
    red_cross = mlines.Line2D([], [], color='red', marker='+', linestyle='None',
                          markersize=10, label='Core & single-copy ortholog')
    green_circle = mlines.Line2D([], [], color='green', marker='.', linestyle='None',
                          markersize=10, label='Accessory')
    magenta_star = mlines.Line2D([], [], color='magenta', marker='*', linestyle='None',
                          markersize=10, label='Unique')

    plt.legend(handles=[blue_x, red_cross, green_circle, magenta_star])

                
    plt.show()