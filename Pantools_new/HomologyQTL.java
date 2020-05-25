/*
 * Detect homologous regions of QTLs in a pangenome.
 * Author: Lily Sakellaridi
 */
package pangenome;

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
//import org.neo4j.graphdb.Transaction;
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
import static pantools.Pantools.disconnect_pangenome;

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
import org.neo4j.graphdb.Result;
//import static pantools.AnnotationLayer;


/**
 *
 * @author sakel001
 */
public class HomologyQTL { // start of class
    public static List<Integer> readGenomes(){
        List<Integer> genomes = new ArrayList<>();
        String filename = "query_genomes.txt";
        String line = null;
        try{
            FileReader fileReader = new FileReader(filename);
            BufferedReader bufferedReader = new BufferedReader(fileReader);
            
            while((line = bufferedReader.readLine()) != null) {
                System.out.println(line);
            }   
            
            bufferedReader.close(); 
        }
            
        catch(FileNotFoundException ex) {
            System.out.println(
                "Unable to open file '" + 
                filename + "'");                
        }
        
        catch(IOException ex) {
            System.out.println(
                "Error reading file '" 
                + filename + "'");                  
            // Or we could just do this: 
            // ex.printStackTrace();
        }
        
        return genomes;
    }
    
    public static List<String> readQTLs() {
        List<String> qtls = new ArrayList<>();
        String filename = "qtls.txt";
        String line = null;
        
        try{
            FileReader fileReader = new FileReader(filename);
            BufferedReader bufferedReader = new BufferedReader(fileReader);
            while((line = bufferedReader.readLine()) != null) {
                System.out.println(line);
            }   
            
            bufferedReader.close(); 
        }
             
        catch(FileNotFoundException ex) {
            System.out.println(
                "Unable to open file '" + 
                filename + "'");                
        }
        
        catch(IOException ex) {
            System.out.println(
                "Error reading file '" 
                + filename + "'");                  
            // Or we could just do this: 
            // ex.printStackTrace();
        }
        
        
        
        return qtls;
    }
    public static void makeHomologyTables(){
        
        List<Integer> genomes = readGenomes();
        List<String> qtls = readQTLs();
        
        
        for (String qtl:qtls){
            
            for (int genome:genomes){
                connect_pangenome();
                String current_query = buildQuery(qtl, genome);
                //System.out.println(current_query);
             
                try(Transaction tx = graphDb.beginTx()){
                    Result res = graphDb.execute(current_query);
                    //String my_res = res.toString();
                    //System.out.print(res);
                    String filename = "Homology table of QTL: " + qtl + " on genome: " + Integer.toString(genome);
                    //System.out.println(filename);
                    writeResultToFile(res, genome, filename);
                 
                tx.success();    
                }
            disconnect_pangenome();      
            }
          
        }
        
    }
    public static void writeResultToFile(Result res, int hom_genome, String filename){
        String i = Integer.toString(hom_genome);
        try {
            FileWriter out = new FileWriter(filename);
            while (res.hasNext()){
                Map<String, Object> row = res.next();
                Object qtl = row.get("q.name");
                Object group = row.get("id(h)");
                Object gene_id = row.get("g.id");
                Object gene_genome = row.get("g.genome");
                Object gene_chr = row.get("g.sequence");
                Object gene_start = row.get("g.begin");
                Object gene_end = row.get("g.end");
                Object gene1_id = row.get("g_h.id");
                Object gene1_genome = row.get("g_h.genome");
                Object gene1_chr = row.get("g_h.sequence");
                Object gene1_start = row.get("g_h.begin");
                Object gene1_end = row.get("g_h.end");
                Object gene2_id = row.get("g" + i + ".id");
                Object gene2_genome = row.get("g" + i + ".genome");
                Object gene2_chr = row.get("g" + i + ".sequence");
                Object gene2_start = row.get("g" + i + ".begin");
                Object gene2_end = row.get("g" + i + ".end");
            
            
                
                System.out.write(qtl);
                System.out.write(",");
                System.out.write(group);
                System.out.write(",");
                System.out.write(gene_id);
                System.out.write(",");
                System.out.write(gene_genome);
                System.out.write(",");
                System.out.write(gene_chr);
                System.out.write(",");
                System.out.write(gene_start);
                System.out.write(",");
                System.out.write(gene_end);
                System.out.write(",");
                System.out.write(gene1_id);
                System.out.write(",");
                System.out.write(gene1_genome);
                System.out.write(",");
                System.out.write(gene1_chr);
                System.out.write(",");
                System.out.write(gene1_start);
                System.out.write(",");
                System.out.write(gene1_end);
                System.out.write(",");
                System.out.write(gene2_id);
                System.out.write(",");
                System.out.write(gene2_genome);
                System.out.write(",");
                System.out.write(gene2_chr);
                System.out.write(",");
                System.out.write(gene2_start);
                System.out.write(",");
                System.out.writeln(gene2_end);
            
            out.close()
            }
            catch(IOException e){
                System.out.println("An error occurred.");
                e.printStackTrace();
            }
            
        
    }
    }
    public static String buildQuery(String qtl, int hom_genome){
        // Builds a Cypher query to make a homology table for a specific QTL-genome combination.
       
        StringBuilder query = new StringBuilder("MATCH (q:QTL) where q.name = \"" + qtl + "\""); // backlashes to wrap the name as string
        query.append("match(m:mRNA {genome: q.genome, sequence: q.sequence})-[:has_homolog]-(h) where (m.begin >= q.begin and m.end <= q.end) ");
        query.append("match(g:gene)-[:codes_for]-(m)");
        
        // Match homologs on same genome, but outside QTL boundaries.
        query.append("optional match(g_h:gene)-[:codes_for]-(m_h:mRNA {genome: q.genome})-[:has_homolog]-(h)");
        query.append("where not (g_h.sequence = q.sequence and g_h.begin >= q.begin and g_h.end <= q.end)");
        
        // Match homologs on the queried genome.
        String i = Integer.toString(hom_genome);
        query.append("optional match(g" + i + ":gene)-[:codes_for]-(m" + i + ":mRNA{genome:" + i + "})-[:has_homolog]-(h)");
        
        // Return result.
        query.append("return q.name, id(h), g.id, g.genome, g.sequence, g.begin, g.end, ");
        query.append("g_h.id, g_h.genome, g_h.sequence, g_h.begin, g_h.end ");    
        query.append(", g" + i + ".id, g" + i + ".genome, g" + i + ".sequence, g" + i + ".begin, g" + i + ".end;");
        
        String query_string = query.toString();
        
        return query_string;
    }
    } // end class HomologyQTL  
