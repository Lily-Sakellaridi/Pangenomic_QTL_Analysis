/*
 * Create a stats summary of QTLs present in a pangenome.
 * Author: Lily Sakellaridi
 */

package pangenome;

// copy imports from EefWP1
import index.kmer;
import pantools.Pantools;
import index.IndexDatabase;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.List;
import java.util.HashMap;
import org.neo4j.graphdb.Label;
import java.util.Scanner;
import java.util.stream.*;
import java.util.stream.Collectors;
import java.util.Map;
import org.apache.commons.lang.ArrayUtils;
import java.util.stream.IntStream;
import index.IndexPointer;
//import org.neo4j.io.fs.FileUtils;
import java.io.InputStreamReader;
import org.neo4j.graphdb.ResourceIterator;
import org.neo4j.graphdb.Transaction;
import org.neo4j.graphdb.GraphDatabaseService;
import org.neo4j.graphdb.Node;
import org.neo4j.graphdb.Relationship;
import org.neo4j.graphdb.RelationshipType;
import org.neo4j.graphdb.ResourceIterator;
import org.neo4j.graphdb.Transaction;
import org.neo4j.graphdb.factory.GraphDatabaseFactory;
import org.neo4j.graphdb.Direction;
import java.math.BigInteger; 
import java.util.concurrent.TimeUnit;

import cern.jet.stat.Gamma;
import org.apache.commons.io.FileUtils;
import static java.lang.Math.toIntExact;
import static pantools.Pantools.PATH_TO_THE_GENOMES_FILE;
import static pantools.Pantools.PATH_TO_THE_REGIONS_FILE;
import static pantools.Pantools.PATH_TO_THE_ANNOTATIONS_FILE;
import static pantools.Pantools.graphDb;
import static pantools.Pantools.homology_group_label;
import static pantools.Pantools.pangenome_label;
import static pantools.Pantools.startTime;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import static pantools.Pantools.write_fasta;
import static pantools.Pantools.total_genomes;
import static pantools.Pantools.adj_total_genomes;
import static pantools.Pantools.target_genome;
import static pantools.Pantools.softcap;
import static pantools.Pantools.softcap_low;
import static pantools.Pantools.WORKING_DIRECTORY;
import static pantools.Pantools.WD_full_path;

import static pantools.Pantools.skip_genomes;
import static pantools.Pantools.PHENOTYPE;
import static pantools.Pantools.GRAPH_DATABASE_PATH;
import static pantools.Pantools.LOG;
import static pantools.Pantools.write_log;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import org.neo4j.graphdb.NotFoundException;
import static pantools.Pantools.phenotype_label;
import static pantools.Pantools.degenerate_label;
import static pantools.Pantools.genome_label;
import static pantools.Pantools.nucleotide_label;
import static pantools.Pantools.gene_label;
import static pantools.Pantools.coding_gene_label;
import static pantools.Pantools.accession_label;
import static pantools.Pantools.mRNA_label;
import static pantools.Pantools.CDS_label;
import static pantools.Pantools.sequence_label;
import static pantools.Pantools.go_label;
import static pantools.Pantools.pfam_label;
import static pantools.Pantools.intron_label;
import static pantools.Pantools.exon_label;
import static pantools.Pantools.tigrfam_label;
import static pantools.Pantools.interpro_label;
import static pantools.Pantools.bgc_label;
import static pantools.Pantools.old_homology_group_label;
import static pantools.Pantools.tRNA_label;
import static pantools.Pantools.rRNA_label;
import static pantools.Pantools.qtl_label; // Lily: add import for QTL label

import static pantools.Pantools.NODE_ID;
import static pantools.Pantools.NODE_ID_long;
import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import static org.neo4j.graphdb.factory.GraphDatabaseSettings.keep_logical_logs;
import static pantools.Pantools.SELECTED_HMGROUPS;
import java.util.HashSet;
import java.util.Set;
import org.neo4j.graphdb.ResourceIterable;
import static org.neo4j.graphdb.factory.GraphDatabaseSettings.keep_logical_logs;
import static pantools.Pantools.GENOME_DATABASE_PATH;
import static pantools.Pantools.GRAPH_DATABASE_PATH;
import static pantools.Pantools.INDEX_DATABASE_PATH;
import static pantools.Pantools.K_SIZE;
import static pantools.Pantools.KMinOne;
import static pantools.Pantools.skip_array;
import static pantools.Pantools.seq_skip_array;
import static pantools.Pantools.seqLayer;
import static pantools.Pantools.OVERWRITE;
import static pantools.Pantools.PATH_TO_THE_PANGENOME_DATABASE;
import pantools.Pantools.RelTypes;
import static pantools.Pantools.pantools_path;
import static pantools.Pantools.current_path;
import static pantools.Pantools.SELECTED_NAME;
import static pantools.Pantools.SELECTED_LABEL;
import static pantools.Pantools.Mode;
import static pantools.Pantools.SELECTED_NAME;
import static pantools.Pantools.NODE_PROPERTY;
import static pantools.Pantools.THREADS;
import static pantools.Pantools.NODE_VALUE;
import static pantools.Pantools.annotation_label;
import static pantools.Pantools.feature_label;
import static pantools.Pantools.feature_label;
import static pantools.Pantools.genomeDb;
import static pantools.Pantools.genome_label;
import static pantools.Pantools.graphDb;
import static pantools.Pantools.indexDb;
import static pantools.Pantools.pangenome_label;
import static pantools.Pantools.startTime;
import sequence.SequenceDatabase;
import sequence.SequenceScanner;
import org.neo4j.graphdb.DatabaseShutdownException;
import index.IndexScanner;
import java.io.BufferedOutputStream;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.Collections;
import java.util.Comparator;
import java.util.Random;
import java.util.TreeSet;
import java.util.concurrent.atomic.AtomicInteger;
import static pantools.Pantools.CONTRAST;
import static pantools.Pantools.INTERSECTION_RATE;
import static pantools.Pantools.MAX_TRANSACTION_SIZE;
import static pantools.Pantools.MCL_INFLATION;
import static pantools.Pantools.MIN_NORMALIZED_SIMILARITY;

import static pantools.Pantools.indexSc;
import static pantools.Pantools.connect_pangenome;

import static pantools.Pantools.phenotype_treshold;
import static pantools.Pantools.geno_pheno_map;
import static pantools.Pantools.new_phenotype_map;
//import static pantools.Pantools.WP3;
import java.text.SimpleDateFormat;  
import java.time.LocalTime;
import java.util.concurrent.ConcurrentHashMap;
import static pangenome.EefWP1.count_nodes_connections;
import static pangenome.EefWP1.include_functional_annotations;
import static pantools.Pantools.go_label;

// Lily; imports for Cypher queries. Doesn work yet.
//import org.neo4j.cypher.CypherParser;
//import org.neo4j.cypher.ExecutionEngine;
//import org.neo4j.cypher.javacompat.ExecutionResult;

/**
 *
 * @author sakel001
 */
public class SummaryQTL {
    public static void printQTLStatistics(String qtlName, int genomeIndex, int sequenceIndex, int qtlSize, int geneCount, int goCount, double lod, double variance, String parentA, String parentB){
        System.out.printf("%s\t%d\t%d\t%d\t%d\t%d\t%f\t%f\t%s\t%s%n", qtlName, genomeIndex, sequenceIndex, qtlSize, geneCount, goCount, lod, variance, parentA, parentB);
        
        
    }
    public static void countGenesPerQTL(){
        System.out.println("Counting genes per QTL...");
        connect_pangenome(); // command: java -jar PATH/pantools.jar qtl_summary -dp PATH_TO_DATABASE
        try (Transaction tx = graphDb.beginTx()){
            Node pangenome_node = graphDb.findNodes(pangenome_label).next();
            total_genomes = (int) pangenome_node.getProperty("num_genomes");
            System.out.printf("The total number of genomes is: %d%n", total_genomes);
            System.out.println("Summary QTL statistics table");
            System.out.println("QTL Name\t Genome\t Sequence (Chromosome or scaffold)\t Number of genes\tNumber of GO terms\tLOD score\tExplained variance\tParent A\tParent B\n");
            //String testCypherQuery = "MATCH (q:QTL) RETURN q";
            //Result result = graphDb.execute(testCypherQuery); does not recognize symbol Result.
            // iterate over genomes
            
            for (int i = 1; i < total_genomes + 1; i++){
                //System.out.println("Genome: " + i);
                ResourceIterator<Node> all_qtl_nodes = graphDb.findNodes(qtl_label, "genome", i);
                
                while (all_qtl_nodes.hasNext()){ 
                    Node qtl_node = all_qtl_nodes.next();
                    //System.out.println(qtl_node);
                    int gene_count = 0; // initialize gene count to zero
                    int go_count = 0;
                    //Iterable<String> all_properties = qtl_node.getPropertyKeys();
                    int qtl_start = (int) qtl_node.getProperty("begin");
                    int qtl_end = (int) qtl_node.getProperty("end");
                    int qtl_size = qtl_end - qtl_start;
                    // the seqid line causes the program to switch to the catch clause and shut down.
                    int qtl_seqid = (int) qtl_node.getProperty("sequence");
                    String qtl_name = (String) qtl_node.getProperty("name");
                    String qtl_parentA = (String) qtl_node.getProperty("parent_a");
                    String qtl_parentB = (String) qtl_node.getProperty("parent_b");
                    double qtl_lod = (double) qtl_node.getProperty("lod_score");
                    double qtl_explained_variance = (double) qtl_node.getProperty("explained_variance");
                    //System.out.printf("QTL %s on genome %d is located on sequence %d%n", qtl_name, i, qtl_seqid);
                    //System.out.printf("It starts at position %d and ends at position %d%n", qtl_start, qtl_end);
                    //System.out.printf("It contains %d genes.", gene_count);
                    // Get all gene nodes that have are in the same genome and sequence as the current QTL
                    ResourceIterator<Node> all_gene_nodes = graphDb.findNodes(gene_label, "genome", i, "sequence", qtl_seqid);
                    while (all_gene_nodes.hasNext()){
                        //System.out.println("Woo hoo");
                        Node gene_node = all_gene_nodes.next();
                        int gene_start = (int) gene_node.getProperty("begin");
                        int gene_end = (int) gene_node.getProperty("end");
                        if (gene_start >= qtl_start && gene_end <= qtl_end){
                            
                            // get the GO terms associated with the current gene
                            // if they are in the go_terms array, ignore them
                            // if they are not, add them
                            
                            //long [] total_go = count_nodes_connections(go_label, RelTypes.has_go);
                            //System.out.println(total_go.length); oh my god it has 41 terms what the hell is in it
                            //System.out.println(total_go);
                            //go_count += total_go[0];
                            /*        
                            if (total_go[0] == 0){
                                include_functional_annotations();
                            } 
                            */
                            
                            gene_count++;
                            
                        }
                    }
                    printQTLStatistics(qtl_name, i, qtl_seqid, qtl_size, gene_count, go_count, qtl_lod, qtl_explained_variance, qtl_parentA, qtl_parentB);
                    
                    
            }
            
        }
            
        }
    
        catch (NotFoundException | ClassCastException no1){
            System.out.println("\nSomething went wrong and stuff.");
        }
       
    }
}
