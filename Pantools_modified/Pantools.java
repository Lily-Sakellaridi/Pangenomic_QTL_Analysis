/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package pantools;

import index.IndexDatabase;
import index.IndexPointer;
import index.IndexScanner;
import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.lang.management.ManagementFactory;
import java.lang.management.MemoryPoolMXBean;
import java.lang.management.MemoryUsage;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.Map;
import java.util.Scanner;
import java.util.concurrent.TimeUnit;
import org.apache.commons.compress.compressors.CompressorInputStream;
import org.apache.commons.compress.compressors.CompressorStreamFactory;
import org.apache.commons.lang.ArrayUtils;
import org.neo4j.graphdb.GraphDatabaseService;
import org.neo4j.graphdb.Label;
import org.neo4j.graphdb.Node;
import org.neo4j.graphdb.Relationship;
import org.neo4j.graphdb.RelationshipType;
import org.neo4j.graphdb.factory.GraphDatabaseFactory;
import static org.neo4j.graphdb.factory.GraphDatabaseSettings.keep_logical_logs;
import org.neo4j.io.fs.FileUtils;
import pangenome.AnnotationLayer;
import pangenome.ProteomeLayer;
import pangenome.GenomeLayer;
import static pangenome.GenomeLayer.append_fwd;
import static pangenome.GenomeLayer.append_rev;
import static pangenome.GenomeLayer.get_outgoing_edge;
import static pangenome.GenomeLayer.locate;
import sequence.SequenceDatabase;
import sequence.SequenceScanner;

//#NEW
import pangenome.EefWP1;
import org.apache.commons.lang3.StringUtils;
import org.neo4j.graphdb.factory.GraphDatabaseSettings;
import pangenome.MLSA;
import pangenome.annotation_pipeline;
//import pangenome.ExecCommand;
import pangenome.EefWP2;
import pangenome.EefWP3;
import pangenome.retired_functions;
//import pangenome.EefWP4;
//import pangenome.NeedlemanWunsch;
//import pangenome.SmithWaterman;

// Lily
import pangenome.SummaryQTL; 
import pangenome.HomologyQTL;

/**
 * Implements the main function and some shared variables and methods. 
 * 
 * @author Siavash Sheikhizadeh, Bioinformatics chairgroup, Wageningen
 * University, the Netherlands
 */
public class Pantools {
    public static String GRAPH_DATABASE_PATH = "/databases/graph.db/";
    public static String INDEX_DATABASE_PATH = "/databases/index.db/";
    public static String GENOME_DATABASE_PATH = "/databases/genome.db/";
    public static String OUTPUT_PATH = "";

    public static GraphDatabaseService graphDb;
    public static IndexDatabase indexDb;
    public static IndexScanner indexSc;    
    public static SequenceDatabase genomeDb;
    public static SequenceScanner genomeSc;    

    public static String WORKING_DIRECTORY;
    public static String WD_full_path;
    public static String PATH_TO_THE_GENOMES_FILE;
    public static String PATH_TO_THE_PROTEOMES_FILE;
    public static String PATH_TO_THE_ANNOTATIONS_FILE;
    public static String PATH_TO_THE_REGIONS_FILE;
    public static String PATH_TO_THE_GENOME_NUMBERS_FILE;
    public static String RAW_ABUNDANCE_FILE = "";
    public static String PATH_TO_THE_FIRST_SRA;
    public static String PATH_TO_THE_SECOND_SRA;
    public static String FEATURE = "gene";
    
    public static boolean CONNECT_ANNOTATIONS = false;
    
    public static double INTERSECTION_RATE = 0.08;
    public static double CONTRAST = 8;
    public static double MCL_INFLATION = 10.8;
    public static int MIN_NORMALIZED_SIMILARITY = 95;

    public static int K_SIZE = -1;
    public static int MAX_ALIGNMENT_LENGTH = 2000;
    public static int GAP_OPEN = -20;
    public static int GAP_EXT = -3;
    public static int ANCHORS_DISTANCE = 10000; // The distance between two anchor nodes
    public static int MAX_TRANSACTION_SIZE = 1000;    //   The number of transactions to be committed in batch
    public static int cores = Runtime.getRuntime().availableProcessors();
    public static long heapSize = Runtime.getRuntime().maxMemory();
    public static boolean DEBUG = false;
    public static boolean SHOW_KMERS;
    public static int THREADS = 1;
    
    public static double MIN_IDENTITY = 0.5;
    public static int NUM_KMER_SAMPLES = 15;
    public static int MAX_NUM_LOCATIONS = 15;
    public static int ALIGNMENT_BOUND = 5;    
    public static int CLIPPING_STRINGENCY = 1; // 0: no-clipping
                                               // 1: low
                                               // 2: medium
                                               // 3: high    
    public static int MIN_HIT_LENGTH = 13;
    public static int MAX_FRAGMENT_LENGTH = 5000;
    public static int SHOULDER = 100;    
    public static int ALIGNMENT_MODE = 2; // 0: all-hits    
                                          // -1: pan-genomic unique_best
                                          // -2: pan-genomic random_best
                                          // -3: pan-genomic all_bests
                                          // 1: genomic unique_best
                                          // 2: genomic random_best
                                          // 3: genomic all_bests
    public static boolean BAMFORMAT = false;
    public static boolean INTERLEAVED = false;
    public static boolean VERYSENSITIVE = false;
    public static boolean SENSITIVE = false;
    public static boolean FAST = true;
    public static boolean VERYFAST = false;
    
    public static Label pangenome_label = Label.label("pangenome");
    public static Label genome_label = Label.label("genome");
    public static Label sequence_label = Label.label("sequence");
    public static Label nucleotide_label = Label.label("nucleotide");
    public static Label degenerate_label = Label.label("degenerate");
    public static Label annotation_label = Label.label("annotation");
    public static Label variation_label = Label.label("variation");
    public static Label gene_label = Label.label("gene");
    public static Label coding_gene_label = Label.label("coding_gene");
    public static Label mRNA_label = Label.label("mRNA");
    public static Label tRNA_label = Label.label("tRNA");
    public static Label rRNA_label = Label.label("rRNA");
    public static Label CDS_label = Label.label("CDS");
    public static Label exon_label = Label.label("exon");
    public static Label intron_label = Label.label("intron");
    public static Label feature_label = Label.label("feature");
    public static Label homology_group_label = Label.label("homology_group");
    public static Label qtl_label = Label.label("QTL"); // Lily: add QTL label
    
    //#NEW all public from Eef
    public static String PATH_TO_THE_PANGENOME_DATABASE;
    public static int total_genomes = 0;
    public static int adj_total_genomes = 0;
    public static int softcap = 0;
    public static int softcap_low = 0;
    public static int KMinOne = 0;
    public static String Mode = "0";
    public static int mismatch = 0;
    public static String NODE_ID;
    public static Long NODE_ID_long;
    public static String pantools_path;
    public static String current_path;
    public static String SELECTED_NAME;
    public static String target_genome;
    public static String skip_genomes;
    public static int [] skip_array;
    public static String [] seq_skip_array;
    public static boolean OVERWRITE = false;
    public static boolean LOG = false;
    public static boolean write_log = false;
    public static boolean UNIQUE_PHENO = false;
    public static String PHENOTYPE;
    public static String SELECTED_LABEL;
    public static String NODE_PROPERTY;
    public static String NODE_VALUE;
    public static Label phenotype_label = Label.label("phenotype");
    public static Label accession_label = Label.label("accession");
    public static String SELECTED_HMGROUPS;
    public static Label go_label = Label.label("GO");
    public static Label interpro_label = Label.label("interpro");
    public static Label bgc_label = Label.label("bcg"); // IMPORTANT! THIS needs to be changed 
    public static Label pfam_label = Label.label("pfam");
    public static Label tigrfam_label = Label.label("tigrfam");
    public static Label old_homology_group_label = Label.label("old_homology_group");
    public static Label low_complexity_label = Label.label("low_complexity");
    public static Label busco_label = Label.label("busco");
    public static Label signalp_label = Label.label("signalp");
    public static Label phobius_label = Label.label("phobius");
    public static Label transmembrane_label = Label.label("transmembrane");
    public static Label progress_label = Label.label("PROGRESS");
    public static Label progress2_label = Label.label("PROGRESS2");
    public static Label rel_node_label = Label.label("rel_node");
    public static Label temp_label = Label.label("temp");
    
    public static HashMap<String, Integer> phenotype_treshold = new HashMap<>();
    public static HashMap<Integer, String> geno_pheno_map;  // genome number is coupled to phenotype  
    public static HashMap<String, int[]> new_phenotype_map = new HashMap<>();
    
    public static enum RelTypes implements RelationshipType {
        FF, FR, RF, RR,
        has, // for pointing to genome and sequence nodes
        starts,
        stops,
        has_homolog, // for pointing to gene nodes from the homology group
        codes_for,// for connecting genes to mRNAs
        is_parent_of,
        contributes_to,// for connecting CDSs to mRNA
        is_similar_to,
        annotates,
        varies, // #new add comma 
        
         //# new eef
        has_phenotype, //#new
        has_go, //#new
        has_tigrfam, //#new
        has_pfam,//#new
        has_interpro,//#new
        is_a,//#new GO
        part_of,//#new GO
        positively_regulates,//#new GO
        negatively_regulates,//#new GO
        regulates,//#new GO
        is_part_of,//#new GO
        has_old_homolog,//#new
        has_busco // #new
    }

    public static char[] sym = new char[]{'A', 'C', 'G', 'T', 'M', 'R', 'W', 'S', 'Y', 'K', 'V', 'H', 'D', 'B', 'N'};
    public static int[] complement = new int[]{3, 2, 1, 0, 9, 8, 6, 7, 5, 4, 13, 12, 11, 10, 14};
    public static int[] binary = new int[256];
   
    public static long startTime;
    public static long phaseTime;
    public static long num_nodes;
    public static int num_degenerates;
    public static long num_edges;
    public static long num_bases;
    public static Node db_node;
    private static String[] label_strings;
    public static Map<String,Label> labels;

    public static GenomeLayer seqLayer;
    public static AnnotationLayer annLayer;
    public static ProteomeLayer proLayer;
    
    /**
     * The main function of PanTools.
     * 
     * @param args The command line arguments
     */
    public static void main(String[] args) {
        int x, i;
        double y;
        File theDir;
        startTime = System.currentTimeMillis();
        if (args.length < 1) {
            print_help_message();
            System.exit(1);
        }
        binary['A'] = 0;
        binary['C'] = 1;
        binary['G'] = 2;
        binary['T'] = 3;
        binary['M'] = 4;
        binary['R'] = 5;
        binary['W'] = 6;
        binary['S'] = 7;
        binary['Y'] = 8;
        binary['K'] = 9;
        binary['V'] = 10;
        binary['H'] = 11;
        binary['D'] = 12;
        binary['B'] = 13;
        binary['N'] = 14; 
        seqLayer = new GenomeLayer();
        annLayer = new AnnotationLayer();
        proLayer = new ProteomeLayer();
        labels = new HashMap<String,Label>();
        label_strings = new String[]{
        "pangenome", "genome","sequence","nucleotide","degenerate",
        "annotation","variation","gene","coding_gene", "mRNA", 
        "tRNA", "rRNA", "CDS", "exon", "intron", "feature", 
        "broken_protein", "homology_group", "low_complexity", "QTL"}; // Lily: add QTL label string        
        for (i = 0; i < label_strings.length; ++i)
            labels.put(label_strings[i], Label.label(label_strings[i]));
        System.out.println("\n------------------------------- PanTools ------------------------------");
        
        //#new eef
        pantools_path = Pantools.class.getProtectionDomain().getCodeSource().getLocation().getPath(); // retrieve location of jar //NEW
        pantools_path = pantools_path.replace("dist/pantools.jar", "");  //NEW
        int loop_counter = 0; //NEW 
        Path currentRelativePath = Paths.get("");  
        current_path = currentRelativePath.toAbsolutePath().toString() + "/";  
        try{
            for (i = 1; i < args.length;){
                loop_counter ++;
                switch (args[i]){
                    case "--kmer-size": case "-ks":
                        x = Integer.parseInt(args[i + 1]);
                        i += 2;
                        if (x >= 6 && x <= 255)
                            K_SIZE = x;
                        else {
                            System.out.println("Choose K in the range [6..255] or do not specify it to be calculated automatically.");
                            System.exit(1);
                        }
                        System.out.println("K_SIZE = " + K_SIZE);
                        break;
                    case "--db-path": case "-dp":
                        WORKING_DIRECTORY = args[i + 1];
                        i += 2;
                        if (!WORKING_DIRECTORY.endsWith("/")){
                            WORKING_DIRECTORY += "/";
                        }
                        System.out.println("WORKING_DIRECTORY = " + WORKING_DIRECTORY);                    
                        PATH_TO_THE_PANGENOME_DATABASE = WORKING_DIRECTORY;
                        create_full_path_working_directory();
                        String [] allowed_commands = new String[]{"bpg","build_pangenome","annotate_genomes","build_panproteome", "bpp" };
                        if (!ArrayUtils.contains(allowed_commands, args[0])){
                            check_database();
                        }                   
                        break;
                    case "--out-path": case "-op":
                        OUTPUT_PATH = args[i + 1];
                        i += 2;
                        theDir = new File(OUTPUT_PATH);
                        if (!theDir.exists()) {
                                System.out.println("OUTPUT_PATH does not exist! Potential outputs will be written in database directory.");
                                OUTPUT_PATH = WORKING_DIRECTORY;
                        } else
                            System.out.println("OUTPUT_PATH = " + OUTPUT_PATH);
                        break;
                    case "--genomes-file": case "-gf":
                        PATH_TO_THE_GENOMES_FILE = args[i + 1];
                        i += 2;
                        theDir = new File(PATH_TO_THE_GENOMES_FILE);
                        if (!theDir.exists()) {
                                System.out.println(PATH_TO_THE_GENOMES_FILE + " does not exist!");
                                System.exit(1);
                        }
                        System.out.println("PATH_TO_THE_GENOMES_FILE = " + PATH_TO_THE_GENOMES_FILE);
                        break;
                    case "--proteomes-file": case "-pf":
                        PATH_TO_THE_PROTEOMES_FILE = args[i + 1];
                        i += 2;
                        theDir = new File(PATH_TO_THE_PROTEOMES_FILE);
                        if (!theDir.exists()) {
                                System.out.println(PATH_TO_THE_PROTEOMES_FILE + " does not exist!");
                                System.exit(1);
                        }
                        System.out.println("PATH_TO_THE_PROTEOMES_FILE = " + PATH_TO_THE_PROTEOMES_FILE);
                        break;
                    case "--annotations-file": case "-af":
                        PATH_TO_THE_ANNOTATIONS_FILE = args[i + 1];
                        i += 2;
                        theDir = new File(PATH_TO_THE_ANNOTATIONS_FILE);
                        if (!theDir.exists()) {
                                System.out.println("\n" + PATH_TO_THE_ANNOTATIONS_FILE + " does not exist! (-af)");
                                System.exit(1);
                        }
                        System.out.println("PATH_TO_THE_ANNOTATIONS_FILE = " + PATH_TO_THE_ANNOTATIONS_FILE);
                        break;
                    case "--connect-annotations": case "-ca":
                        CONNECT_ANNOTATIONS = true;
                        i += 1;
                        System.out.println("CONNECT_ANNOTATIONS = true");
                        break;
                    case "--regions-file": case "-rf":
                        PATH_TO_THE_REGIONS_FILE = args[i + 1];
                        i += 2;
                        theDir = new File(PATH_TO_THE_REGIONS_FILE);
                        if (!theDir.exists()) {
                                System.out.println(PATH_TO_THE_REGIONS_FILE + " does not exist!");
                                System.exit(1);
                        }
                        System.out.println("PATH_TO_THE_REGIONS_FILE = " + PATH_TO_THE_REGIONS_FILE);
                        break;
                    case "--genome-numbers": case "-gn":
                        PATH_TO_THE_GENOME_NUMBERS_FILE = args[i + 1];
                        i += 2;
                        theDir = new File(PATH_TO_THE_GENOME_NUMBERS_FILE);
                        if (!theDir.exists()) {
                                System.out.println(PATH_TO_THE_GENOME_NUMBERS_FILE + " does not exist!");
                                System.exit(1);
                        }
                        System.out.println("PATH_TO_THE_GENOME_NUMBERS_FILE = " + PATH_TO_THE_GENOME_NUMBERS_FILE);
                        break;
                    case "--intersection-rate": case "-ir": 
                        y = Double.parseDouble(args[i + 1]);
                        i += 2;
                        if (y >= 0.001 && y <= 0.1)
                            INTERSECTION_RATE = y;
                        else {
                            System.out.println("Choose INTERSECTION_RATE in the range [0.001..0.1] or do not specify it to use the default value.");
                            System.exit(1);
                        }
                        System.out.println("INTERSECTION_RATE = " + INTERSECTION_RATE);
                        break;
                    case "--similarity-threshold": case "-st": 
                        x = Integer.parseInt(args[i + 1]);
                        i += 2;
                        if (x > 0 && x < 100)
                            MIN_NORMALIZED_SIMILARITY = x;
                        else {
                            System.out.println("Choose MIN_NORMALIZED_SIMILARITY in the range ]0..100[ or do not specify it to use the default value.");
                            System.exit(1);
                        }
                        System.out.println("MIN_NORMALIZED_SIMILARITY = " + MIN_NORMALIZED_SIMILARITY);
                        break;
                    case "--mcl-inflation": case "-mi": 
                        y = Double.parseDouble(args[i + 1]);
                        i += 2;
                        if (y > 1 && y < 19)
                            MCL_INFLATION = y;
                        else {
                            System.out.println("Choose MCL_INFLATION in the range ]1.0..19.0[ or do not specify it to use the default value.");
                            System.exit(1);
                        }
                        System.out.println("MCL_INFLATION = " + MCL_INFLATION);
                        break;
                    case "--contrast": case "-ct": 
                        y = Double.parseDouble(args[i + 1]);
                        i += 2;
                        if (y > 0 && y < 10)
                            CONTRAST = y;
                        else {
                            System.out.println("Choose CONTRAST in the range [1..9] or do not specify it to use the default value.");
                            System.exit(1);
                        }
                        System.out.println("CONTRAST = " + CONTRAST);
                        break;
                    case "--relaxation": case "-rn": 
                        x = Integer.parseInt(args[i + 1]);
                        i += 2;
                        if (x >= 1 && x <= 8){
                            INTERSECTION_RATE = new double[] {0, 0.08, 0.07, 0.06, 0.05, 0.04, 0.03, 0.02, 0.01}[x];
                            MIN_NORMALIZED_SIMILARITY = new int[]   {0,95, 85, 75, 65, 55, 45, 35, 25}[x];
                            MCL_INFLATION = new double[]{0, 10.8, 9.6, 8.4, 7.2, 6.0, 4.8, 3.6, 2.4}[x];
                            CONTRAST = new double[] {0,8, 7, 6, 5, 4, 3, 2, 1 }[x];
                        }
                        else {
                            System.out.println("Choose RELAXATION in the range [1..8] or do not specify it to use the default values.");
                            System.exit(1);
                        }
                        System.out.println("INTERSECTION = " + INTERSECTION_RATE);
                        System.out.println("MIN_NORMALIZED_SIMILARITY = " + MIN_NORMALIZED_SIMILARITY);
                        System.out.println("MCL_INFLATION = " + MCL_INFLATION);
                        System.out.println("CONTRAST = " + CONTRAST);
                        break;
                     
                    case "--threads_number": case "-tn":
                        x = Integer.parseInt(args[i + 1]);
                        i += 2;
                        if (x < cores)
                            THREADS = x;
                        else {
                            System.out.println("The maximum number of threads on this machine = " + cores + ".");
                            THREADS = cores;
                        }
                        System.out.println("THREADS = " + THREADS);
                        break;
                    case "--gap-open": case "-go":
                        x = Integer.parseInt(args[i + 1]);
                        i += 2;
                        if (x >= -50 && x <= -1)
                            GAP_OPEN = x;
                        else {
                            System.out.println("Choose GAP_OPEN in the range [-50..-1] or do not specify it to use the default value.");
                            System.exit(1);
                        }
                        System.out.println("GAP_OPEN = " + GAP_OPEN);
                        break;
                    case "--gap-extention": case "-ge":
                        x = Integer.parseInt(args[i + 1]);
                        i += 2;
                        if (x >= -5 && x <= -1)
                            GAP_EXT = x;
                        else {
                            System.out.println("Choose GAP_EXT in the range [-5..-1] or do not specify it to use the default value.");
                            System.exit(1);
                        }
                        System.out.println("GAP_EXT = " + GAP_EXT);
                        break;
                    case "--feature_type": case "-ft":
                        if (labels.containsKey(args[i + 1])){
                            FEATURE = args[i + 1];
                            i += 2;
                        } else {
                            System.out.println(args[i + 1] + " is not a unknown feature.");
                            System.exit(1);
                        }
                        System.out.println("FEATURE = " + FEATURE);
                        break;
                    case "--first_sra": case "-1":
                        PATH_TO_THE_FIRST_SRA = args[i + 1];
                        i += 2;
                        theDir = new File(PATH_TO_THE_FIRST_SRA);
                        if (!theDir.exists()) {
                                System.out.println(PATH_TO_THE_FIRST_SRA + " does not exist!");
                                System.exit(1);
                        }
                        System.out.println("PATH_TO_THE_FIRST_SRA = " + PATH_TO_THE_FIRST_SRA);
                        break;
                    case "--second_sra": case "-2":
                        PATH_TO_THE_SECOND_SRA = args[i + 1];
                        i += 2;
                        theDir = new File(PATH_TO_THE_SECOND_SRA);
                        if (!theDir.exists()) {
                                System.out.println(PATH_TO_THE_SECOND_SRA + " does not exist!");
                                System.exit(1);
                        }
                        System.out.println("PATH_TO_THE_SECOND_SRA = " + PATH_TO_THE_SECOND_SRA);
                        break;
                    case "--clip_strigency": case "-cs":
                        CLIPPING_STRINGENCY = Integer.parseInt(args[i + 1]);
                        i += 2;
                        switch (CLIPPING_STRINGENCY){
                            case 0:
                            System.out.println("CLIPPING_STRINGENCY = " + CLIPPING_STRINGENCY + " : no-clipping");
                                break;
                            case 1:
                            System.out.println("CLIPPING_STRINGENCY = " + CLIPPING_STRINGENCY + " : low");
                                break;
                            case 2:
                            System.out.println("CLIPPING_STRINGENCY = " + CLIPPING_STRINGENCY + " : medium");
                                break;
                            case 3:
                            System.out.println("CLIPPING_STRINGENCY = " + CLIPPING_STRINGENCY + " : high");
                                break;
                            default:
                                System.out.println("Choose CLIPPING_STRINGENCY 0, 1, 2, or 3, or do not specify it to use the default value of 0.");
                                System.exit(1);
                        }
                        break;
                    case "--very-sensitive": case "-vs":
                        VERYSENSITIVE = true;
                        MIN_IDENTITY = 0.5;
                        NUM_KMER_SAMPLES = 30;
                        MAX_NUM_LOCATIONS = 30;
                        ALIGNMENT_BOUND = 12;    
                        CLIPPING_STRINGENCY = 3;
                        i += 1;
                        System.out.println("MIN_IDENTITY = 0.5");
                        System.out.println("NUM_KMER_SAMPLES = 30");
                        System.out.println("MAX_NUM_LOCATIONS = 30");
                        System.out.println("ALIGNMENT_BOUND = 12");
                        System.out.println("CLIPPING_STRINGENCY = 3");
                        break;
                    case "--sensitive": case "-sv":
                        SENSITIVE = true;
                        MIN_IDENTITY = 0.6;
                        NUM_KMER_SAMPLES = 23;
                        MAX_NUM_LOCATIONS = 23;
                        ALIGNMENT_BOUND = 9;    
                        CLIPPING_STRINGENCY = 2;
                        i += 1;
                        System.out.println("MIN_IDENTITY = 0.6");
                        System.out.println("NUM_KMER_SAMPLES = 23");
                        System.out.println("MAX_NUM_LOCATIONS = 23");
                        System.out.println("ALIGNMENT_BOUND = 9");
                        System.out.println("CLIPPING_STRINGENCY = 2");
                        break;
                    case "--fast": case "-f":
                        FAST = true;
                        MIN_IDENTITY = 0.7;
                        NUM_KMER_SAMPLES = 15;
                        MAX_NUM_LOCATIONS = 15;
                        ALIGNMENT_BOUND = 6;    
                        CLIPPING_STRINGENCY = 1;
                        i += 1;
                        System.out.println("MIN_IDENTITY = 0.7");
                        System.out.println("NUM_KMER_SAMPLES = 15");
                        System.out.println("MAX_NUM_LOCATIONS = 15");
                        System.out.println("ALIGNMENT_BOUND = 6");
                        System.out.println("CLIPPING_STRINGENCY = 1");
                        break;
                    case "--very-fast": case "-vf":
                        VERYFAST = true;
                        MIN_IDENTITY = 0.8;
                        NUM_KMER_SAMPLES = 7;
                        MAX_NUM_LOCATIONS = 7;
                        ALIGNMENT_BOUND = 3;    
                        CLIPPING_STRINGENCY = 0;
                        i += 1;
                        System.out.println("MIN_IDENTITY = 0.8");
                        System.out.println("NUM_KMER_SAMPLES = 7");
                        System.out.println("MAX_NUM_LOCATIONS = 7");
                        System.out.println("ALIGNMENT_BOUND = 3");
                        System.out.println("CLIPPING_STRINGENCY = 0");
                        break;
                    case "--interleaved": case "-il":
                        INTERLEAVED = true;
                        i += 1;
                        System.out.println("INTERLEAVED = true");
                        break;
                    case "--bam-format": case "-bf":
                        BAMFORMAT = true;
                        i += 1;
                        System.out.println("BAMFORMAT = true");
                        break;
                    case "--min_mapping-identity": case "-mmi":
                        y = Double.parseDouble(args[i + 1]);
                        i += 2;
                        if (y >= 0 && y < 1)
                           MIN_IDENTITY = y;
                        else {
                            System.out.println("Choose MIN_IDENTITY in the range [0..1[ or do not specify it to use the default value.");
                            System.exit(1);
                        }
                        System.out.println("MIN_IDENTITY = " + MIN_IDENTITY);
                        break;
                    case "--num-kmer-samples": case "-nks":
                        x = Integer.parseInt(args[i + 1]);
                        i += 2;
                        if (x >= 1)
                           NUM_KMER_SAMPLES = x;
                        else {
                            System.out.println("Choose a non-zero NUM_KMER_SAMPLES or do not specify it to use the default value.");
                            System.exit(1);
                        }
                        System.out.println("NUM_KMER_SAMPLES = " + NUM_KMER_SAMPLES);
                        break;
                    case "--max-alignment-length": case "-mal":
                        x = Integer.parseInt(args[i + 1]);
                        i += 2;
                        if (x >= 50 && x<=5000)
                           MAX_ALIGNMENT_LENGTH = x;
                        else {
                            System.out.println("Choose MAX_ALIGNMENT_LENGTH in the range [50..5000] or do not specify it to use the default value.");
                            System.exit(1);
                        }
                        System.out.println("MAX_ALIGNMENT_LENGTH = " + MAX_ALIGNMENT_LENGTH);
                        break;
                    case "--max-fragment-length": case "-mfl":
                        x = Integer.parseInt(args[i + 1]);
                        i += 2;
                        if (x >= 50 && x<=5000)
                           MAX_FRAGMENT_LENGTH = x;
                        else {
                            System.out.println("Choose MAX_FRAGMENT_LENGTH in the range [50..5000] or do not specify it to use the default value.");
                            System.exit(1);
                        }
                        System.out.println("MAX_FRAGMENT_LENGTH = " + MAX_FRAGMENT_LENGTH);
                        break;
                    case "--min-hit_length": case "-mhl":
                        x = Integer.parseInt(args[i + 1]);
                        i += 2;
                        if (x >= 10 && x <= 100)
                           MIN_HIT_LENGTH = x;
                        else {
                            System.out.println("Choose MIN_HIT_LENGTH in the range [10..100] or do not specify it to use the default value.");
                            System.exit(1);
                        }
                        System.out.println("MIN_HIT_LENGTH = " + MIN_HIT_LENGTH);
                        break;
                    case "--alignment_bound": case "-ab":
                        x = Integer.parseInt(args[i + 1]);
                        i += 2;
                        if (x >= 1 && x <= 100)
                           ALIGNMENT_BOUND = x;
                        else {
                            System.out.println("Choose ALIGNMENT_BOUND in the range [1..100] or do not specify it to use the default value.");
                            System.exit(1);
                        }
                        System.out.println("ALIGNMENT_BOUND = " + ALIGNMENT_BOUND);
                        break;
                    case "--max-num-locations": case "-mnl":
                        x = Integer.parseInt(args[i + 1]);
                        i += 2;
                        if (x >= 1 && x <= 100)
                           MAX_NUM_LOCATIONS = x;
                        else {
                            System.out.println("Choose MAX_NUM_LOCATIONS in the range [1..100] or do not specify it to use the default value.");
                            System.exit(1);
                        }
                        System.out.println("MAX_NUM_LOCATIONS = " + MAX_NUM_LOCATIONS);
                        break;
                    case "--alignment-mode": case "-am":
                        ALIGNMENT_MODE = Integer.parseInt(args[i + 1]);
                        i += 2;
                        switch (ALIGNMENT_MODE){
                            case 0:
                                System.out.println("ALIGNMENT_MODE = " + ALIGNMENT_MODE + " : all-hits");
                            break;
                            case -1:
                                System.out.println("ALIGNMENT_MODE = " + ALIGNMENT_MODE + " : unique pangenomic-best");
                            break;
                            case -2:
                                System.out.println("ALIGNMENT_MODE = " + ALIGNMENT_MODE + " : random pangenomic-best");
                            break;
                            case -3:
                                System.out.println("ALIGNMENT_MODE = " + ALIGNMENT_MODE + " : all pangenomic-bests");
                            break;
                            case 1:
                                System.out.println("ALIGNMENT_MODE = " + ALIGNMENT_MODE + " : unique genomic-best");
                            break;
                            case 2:
                                System.out.println("ALIGNMENT_MODE = " + ALIGNMENT_MODE + " : random genomic-best");
                            break;
                            case 3:
                                System.out.println("ALIGNMENT_MODE = " + ALIGNMENT_MODE + " : all genomic-bests");
                            break;
                            default:    
                                System.out.println("Choose ALIGNMENT_MODE in range [-3..3] or leave it to use the default value.");
                                System.exit(1);
                        }
                        break;
                    case "--raw-abundance-file": case "-raf":
                        RAW_ABUNDANCE_FILE = args[i + 1];
                        i += 2;
                        theDir = new File(RAW_ABUNDANCE_FILE);
                        if (!theDir.exists()) {
                                System.out.println(RAW_ABUNDANCE_FILE + " does not exist!");
                                System.exit(1);
                        }
                        System.out.println("RAW_ABUNDANCE_FILE = " + RAW_ABUNDANCE_FILE);
                        break;
                        
                    
                        
                    //#NEW Eef
                    case "--phenotype": case "-p": case "-phenotype": case "-pheno":
                        PHENOTYPE = args[i + 1];
                        System.out.println("SELECTED PHENOTYPE = '" + PHENOTYPE + "'");
                        i += 2;
                        break;
                    case "--property": case "-property":
                        NODE_PROPERTY = args[i + 1];
                        i += 2;
                        break;
                    case "--value": case "-value":
                        NODE_VALUE = args[i + 1];
                        i += 2;
                        break;
                    case "--label": case "-label":
                        SELECTED_LABEL = args[i + 1];
                        i += 2;
                        break;
                    case "-node": case "--node":          
                        NODE_ID = args[i + 1];   
                        System.out.println("NODE = " + NODE_ID);
                        if(NODE_ID.matches("^[0-9]*$") && NODE_ID.length() > 0){ // if string ONLY contains numbers 
                             NODE_ID_long = Long.parseLong(NODE_ID);
                        }
                        i += 2;
                        break;
                    case "-name": case "--name":
                        SELECTED_NAME = args[i + 1];
                        i += 2;
                        break;
                    case "-mode": case "--mode":
                        Mode = args[i + 1];
                        Mode = Mode.toUpperCase();
                        System.out.println("MODE = " + Mode);
                        i += 2;
                        break;
                    case "-softcap_high": case "--softcap_high":
                        String softcap_high_str = args[i + 1];
                        softcap = Integer.parseInt(softcap_high_str);
                        i += 2;
                        break;
                    case "-softcap_low": case "--softcap_low":
                        String softcap_low_str = args[i + 1];
                        softcap_low = Integer.parseInt(softcap_low_str);
                        i += 2;
                        break;    
                    case "--reference": case "-reference": case "-ref": case "--ref":
                        target_genome  = args[i + 1];
                        i += 2;
                        break;
                    case "--homology": case "-homology": case "-hm": case "--hm":
                        SELECTED_HMGROUPS = args[i + 1];     
                        System.out.println("Homology groups = " + SELECTED_HMGROUPS);
                        i += 2;
                        break;
                    case "--skip":case "-skip":
                        skip_genomes = args[i + 1];
                        i += 2;
                        break;
                    case "--log":case "-log":
                        System.out.println("Logging extra output = true");
                        LOG = true;
                        i++;
                        break;
                    case "--eef": case "-eef":
                        print_eeftools(); //NEW
                        i++;
                        break;
                    case "--overwrite":case "-overwrite":
                        OVERWRITE = true;
                        i++;
                        break; 
                    case "--help": case "-h":
                        print_help_message();
                        System.exit(1);
                        break;
                }  
                
                if (loop_counter > 50){
                      System.out.println("\nUnable to parse the given arguments: '" + args[i] + "'\n");
                      
                      System.exit(1);
                }
              
            }
        } catch (NumberFormatException ex){
            System.out.println("The given number is not in the correct format!");
            System.exit(1);
        }
        if (OUTPUT_PATH.equals("")){
            //System.out.println("OUTPUT_PATH has not been specified. Potential outputs will be written in " + WORKING_DIRECTORY + ".");
            OUTPUT_PATH = WORKING_DIRECTORY;
        }
        switch (args[0]) {
            
             case "annotate_genomes":
                annotation_pipeline.annotate_genomes();
                break;
            
            // Lily
            case "qtl_summary":
                SummaryQTL.countGenesPerQTL();
                break;
            case "detect_homologous_qtls":
                HomologyQTL.makeHomologyTables();
                break;
            
            //WP1
            case "compare_grouping":  
                EefWP1.compare_grouping();
                break;
            case "change_active_hmgroups": case "change_active_homology_groups":
                EefWP1.change_active_homology_groups();
                break;
             case "ancestral_reconstruction": 
                EefWP1.ancestral_reconstruction();
                break;
            case "homology_group_overview": case "homology_groups_overview": case "hmgroup_overview":
                EefWP1.homology_group_overview();
                break;
            case "function_of_mRNAs": case "function_of_mrnas":
                EefWP1.function_of_mRNAs();
                break;
            case "functional_annotation_overview": case "FA_overview": case "fa_overview":
                EefWP1.functional_annotation_overview();
                break;
            case "add_phenotype":
                EefWP1.add_phenotype();
                break;
            case "remove_nodes":
                EefWP1.remove_nodes();
                break;
            case "remove_node_property":
                EefWP1.remove_node_property();
                break;
            case "gene_classification":
                EefWP1.gene_classification();
                break;
            case "functional_classification":
                EefWP1.functional_classification();
                break;
            case "kmer_frequency":
                EefWP1.kmer_frequency();
                break;
            case "go_enrichment": case "GO_enrichment":
                EefWP1.go_enrichment();
                break;
            case "show_hmgroup": case "show_homology_group": case "function_of_hmgroup": case "function_of_homology_group":
                EefWP1.show_homology_group();
                break;         
            case "hmgroup_differences_function":
                EefWP1.hmgroup_differences_function();
                break;
            case "create_input_files": // creates genome and annotation input files for the construction steps 
                EefWP1.create_input_files();
                break;
            case "add_label":
                EefWP1.add_label();
                break;
            case "score_phylogeny":
                EefWP1.score_phylogeny();
                break; 
            case "add_property":
                EefWP1.add_property();
                break;
            case "remove_label":
                EefWP1.remove_label();
                break;
            case "pangenome_size_genes": case "pangenome_size_gene" : case "pangenome_size":
                EefWP1.pangenome_size();
                break;  
            case "pangenome_size_kmer": case "pangenome_size_kmers":
                EefWP1.pangenome_size_kmer();
                break;  
            case "add_functional_annotations":
                EefWP1.add_interpro_results();
                break;
             case "report_genes":
                EefWP1.report_genes();
                break;  
             case "signalp":
                EefWP1.signalp();
                break;  
            case "phobius":
                EefWP1.phobius();
                break;  
            case "kmer_classification":
                EefWP1.kmer_classification(); 
                break;
            case "find_genome_number":
                EefWP1.find_genome_number_standalone();
                break;
            case "compare_bgc":
                EefWP1.compare_bgc();
                break;
            case "add_antismash_results":
                EefWP1.add_antismash_results();
                break;
            case "add_homology_groups": case "add_hmgroups":
                EefWP1.add_homology_groups();
                break;
            case "phenotype_overview":
                EefWP1.show_phenotypes(true);
                break;
            case "calculate_ani":
                EefWP1.calculate_ani();
                break;
            case "create_tree_templates": case "create_tree_template":
                EefWP1.create_tree_templates();
                break;
            case "change_genome_paths":
                EefWP1.change_genome_paths();
                break;
            case "locate_genes":
                EefWP1.locate_genes();
                break;
            case "show_go":
                EefWP1.show_go();
                break;
            case "genome_overview":
                EefWP1.show_genomes();
                break;
            case "show_annotations": case "annotation_overview":
                EefWP1.show_annotations();
                break;
            case "find_genes_by_annotation":
                EefWP1.find_genes_by_annotation();
                break; 
            case "msa_of_hmgroups":
                System.out.println("\nThe function 'msa_of_hmgroups' does not exist. Did u mean 'msa_of_hmgroup'?\n");
                break;
            case "msa_of_hmgroup":
                EefWP1.msa_of_hmgroup(true, false, false);
                //EefWP1.msa_of_hmgroup(false, true, true);
                break;
            case "msa_of_hmgroup_no_trimming":
                EefWP1.msa_of_hmgroup_no_trimming();
                break;
             case "most_similar_kmers":  
                EefWP1.most_similar_kmers();
                break;
            case "most_similar_genes":
                EefWP1.most_similar_genes();
                break;
            case "show_node":
                EefWP1.show_node(true);
                break; 
            //case "most_similar_core_snps_kmers":
            //    EefWP1.most_similar_core_snps_kmers();
            //   break;
            case "most_similar_sco_snps":
                EefWP1.most_similar_sco_snps();
                break;
            case "species_tree":
                EefWP1.species_tree();
                break;
            case "basic_statistics": case "statistics":
                EefWP1.basic_statistics();
                break;
            case "move_hmgroups": case "move_homology_groups":
                EefWP1.move_hmgroups(true);
                break;
            case "remove_hmgroup": case "remove_homology_group":
                System.out.println("\n'" + args[0] + "' does not exist! Did u mean 'remove_hmgroups' or 'remove_homology_groups'?\n");
                break;
            case "remove_hmgroups": case "remove_homology_groups":
                EefWP1.remove_hmgroups();
                break;
            case "show_genes":
                EefWP1.show_genes();
                break;
            case "query":
                EefWP1.query();
                break;
            case "blast":
                EefWP1.blast();
                break;
            case "find_genes_in_region":
                EefWP1.find_genes_in_region();
                break;
            case "core_accessory_range":
                EefWP1.core_accessory_range();
                break;
            case "compare_go":
                EefWP1.compare_go();
                break;
            case "msa_of_regions":
                EefWP1.msa_of_regions();
                break;       
            case "mlsa_find_genes": case "MLSA_find_genes":
                MLSA.mlsa_find_genes();
                break;
            case "mlsa_concatenate": case "MLSA_concatenate":
                MLSA.mlsa_concatenate();
                break;
            case "mlsa_check_sequences": case "MLSA_check_sequences":   
                MLSA.mlsa_check_sequences();
                break;
            case "MLSA": case "mlsa":  
                MLSA.run_MLSA();
                break;
             
            //WP2
            case "add_accessions":
                EefWP2.add_accessions();
                break;
            case "find_homologous_regions":
                EefWP2.find_homologous_regions();
                break;
            case "compare_mums":
                EefWP2.compare_mums();
                break;     
            case "run_mummer":
                EefWP2.run_mummers();
                break;
            case "construct_mums":
                EefWP2.construct_mums();
                break;
                 
            //WP3
            case "busco_protein":
                EefWP3.busco_protein();
                break;
            case "busco_genome":
                EefWP3.busco_genome();
                break; 
            case "busco_overview":
                EefWP3.busco_overview();
                break; 
            // case "check_uniprot":
            //    EefWP1.check_uniprot(); // not finished
            //    break;
            case "find_kmer":
                 EefWP1.find_kmer();
                 break;
            case "TEST":
                 EefWP3.test_database();
                 break;
                             
            // for development only
            case "validate_hmgroups":
                 proLayer.validate_hmgroups();
                 break;       
                 
            case "remove_coordinates":
                seqLayer.remove_coordinates();
                break;     
                 
            // retired functions. These may still work or maybe not 
            case "hmgroup_differences_sequence":  // use msa_of_hmgroup instead of this function 
                retired_functions.hmgroup_differences_sequence();
                break;          
            case "one_vs_all":  
                retired_functions.one_vs_all();
                break;
            case "show_phenotype_row":
                EefWP1.show_phenotype_row();
                break;
            case "show_gene_addresses":
                EefWP1.show_gene_addresses();
                break;   
            //case "calc_jaccard_dist_shared_kmers":
            //    EefWP1.calc_jaccard_dist_shared_kmers();
            //    break;    
            // end of retired functions
            
            case "print_nucleotide_nodes":
                 EefWP1.print_nucleotide_nodes();
                break;
            case "bpg": case "build_pangenome":
                seqLayer.initialize_pangenome(false);
                break;
            case "build_testgenome":
                seqLayer.initialize_pangenome(true);
                break;
             case "localize_nodes":
                seqLayer.localize_test();
                break;
            case "bpp": case "build_panproteome":
                proLayer.initialize_panproteome();
                break;
            case "ag": case "add_genomes":
                seqLayer.add_genomes();
                break;
            case "aa": case "add_annotations":
                annLayer.add_annotations();
                break;
            case "ra": case "remove_annotations":
                annLayer.remove_annotaions();
                break;
            case "g": case "group":
                //INTERSECTION_RATE = 0.06;
                //MIN_NORMALIZED_SIMILARITY = 75;
                MCL_INFLATION = 1.5;  
                //CONTRAST  = 6;
                proLayer.group();
                
                // NEW Eef

                break;
           case "find_best_grouping":
                proLayer.run_all_group_settings();
                EefWP1.compare_busco_to_grouping();
                
                break;
            case "rh": case "remove_homologies":
                proLayer.remove_homology_groups();
                break;
            case "rf": case "retrieve_features":
                annLayer.retrieve_feature();
                break;
            case "rr": case "retrieve_regions":
                seqLayer.retrieve_regions();
                break;
            case "rg": case "retrieve_genomes":
                seqLayer.retrieve_genomes();
                break;
            case "rs": case "retrieve_synteny":
                seqLayer.retrieve_synteny(args[2]);
                break;
            case "m": case "map":
                seqLayer.map_reads();
                break;
            case "h": case "help": case "-h": case "--help":
                print_help_message();
                System.exit(1);
                break;
            case "v": case "version":
                System.out.println("PanTools version 2.0\nNeo4j community edition 3.5.1\nJDK 1.8");
                System.exit(1);
            default:
                System.out.println(args[0] + " is not a valid PanTools command, type 'pantools.jar h [or help]' to see the mannual.");
                System.exit(1);
        }
        System.out.println("Total time : " + (System.currentTimeMillis() - startTime) / 1000 + "." + (System.currentTimeMillis() - startTime) % 1000 + " seconds");
        print_peak_memory();
        System.out.println("-----------------------------------------------------------------------");
    }

    /**
     * Connects to genome, index and graph databases of the pan-genome.
     */
    public static void connect_pangenome(){
        if (WORKING_DIRECTORY == null){
            System.out.println("No database was provided (-dp)\n");
            System.exit(1);
        }
        Scanner s;
        String str;
        if (new File(WORKING_DIRECTORY).exists()){
            if (! new File(WORKING_DIRECTORY + GRAPH_DATABASE_PATH).exists()) {
                System.out.println("No graph database found at " + WORKING_DIRECTORY);
                System.exit(1);
            }
            if (! new File(WORKING_DIRECTORY + INDEX_DATABASE_PATH).exists()) {
                System.out.println("No index database found at " + WORKING_DIRECTORY);
                System.exit(1);
            }
            graphDb = new GraphDatabaseFactory().newEmbeddedDatabaseBuilder(
                 new File(WORKING_DIRECTORY + GRAPH_DATABASE_PATH))
                .setConfig(keep_logical_logs, "4 files")
                //.setConfig( GraphDatabaseSettings.pagecache_memory, "5012M" ) 
                // .setConfig( GraphDatabaseSettings.pagecache_memory, "512M" ) 
                .newGraphDatabase();  
            registerShutdownHook(graphDb);
            indexDb = new IndexDatabase(WORKING_DIRECTORY + INDEX_DATABASE_PATH, "sorted");
            if (! new File(WORKING_DIRECTORY + GENOME_DATABASE_PATH).exists()){
                s = new Scanner(System.in);
                System.out.println("No genome database found at " + WORKING_DIRECTORY);
                System.out.println("Do you want to reconstruct it from the graph database [y/n]? ");
                str = s.nextLine().toLowerCase();
                while (!str.equals("y") && !str.equals("n")){
                System.out.println("Do you want to reconstruct it from the graph database [y/n]? ");
                     str = s.nextLine().toLowerCase();
                }
                if (str.equals("y")){
                    rebuild_genome_database();
                } else {
                    System.out.println("Exiting the program...");
                    System.exit(1);  
                }
            } else
                genomeDb = new SequenceDatabase(WORKING_DIRECTORY + GENOME_DATABASE_PATH);
        } else {
            System.out.println("No pangenome found at " + WORKING_DIRECTORY);
            System.exit(1);
        }
    }    
    
    /**
     * Disconnects genome, index and graph databases of the pan-genome.
     */
    public static void disconnect_pangenome(){
        graphDb.shutdown();
        genomeDb.close();
        indexDb.close();
        File directory = new File(WORKING_DIRECTORY + GRAPH_DATABASE_PATH);
        for (File f : directory.listFiles()) {
            if (f.getName().startsWith("neostore.transaction.db.")) {
                f.delete();
            }
        }
    }    
    
    /**
     * Creates and connects to genome, index and graph databases of the pan-genome.
     */
    public static void create_pangenome_database(){
        File theDir;
        Scanner s;
        String str;
        s = new Scanner(System.in);
        if (WORKING_DIRECTORY == null){
            System.out.println("WORKING_DIRECTORY is empty.");
            System.exit(1);
        }
        theDir = new File(WORKING_DIRECTORY);
        if (theDir.exists()) {
            System.out.println("A pangenome database already exists at " + WORKING_DIRECTORY + GRAPH_DATABASE_PATH + ".");
            System.out.println("Do you want to connect to it? otherwisw it would be removed [y/n]? ");
            str = s.nextLine().toLowerCase();
            while (!str.equals("y") && !str.equals("n")){
                 System.out.println("Do you want to connect to it? otherwisw it would be removed [y/n]? ");
                 str = s.nextLine().toLowerCase();
            }
            if (str.equals("y")) {
                connect_pangenome();
                return;
            } else {
                try {
                    FileUtils.deleteRecursively(new File(WORKING_DIRECTORY + GRAPH_DATABASE_PATH));
                } catch (IOException ioe) {
                    System.out.println("Failed to delete the graph database");
                    System.exit(1);  
                }
            }
        }
        try {
            theDir.mkdir();
        } catch (SecurityException se) {
            System.out.println("Failed to create directory " + WORKING_DIRECTORY);
            System.exit(1);
        }
        graphDb = new GraphDatabaseFactory().newEmbeddedDatabaseBuilder(new File(WORKING_DIRECTORY + GRAPH_DATABASE_PATH))
            .setConfig(keep_logical_logs, "4 files").newGraphDatabase();  
        registerShutdownHook(graphDb);
        genomeDb = new SequenceDatabase(WORKING_DIRECTORY + GENOME_DATABASE_PATH, PATH_TO_THE_GENOMES_FILE);
        indexDb = new IndexDatabase(WORKING_DIRECTORY + INDEX_DATABASE_PATH, PATH_TO_THE_GENOMES_FILE, genomeDb, K_SIZE);
        EefWP1.create_directory_in_DB_dir("/databases/genome.db/genomes");
    }

    /**
     * Creates and connects to graph databases of the pan-genome.
     */
    public static void create_panproteome_database(){
        File theDir;
        Scanner s;
        String str;
        s = new Scanner(System.in);
        if (WORKING_DIRECTORY == null){
            System.out.println("WORKING_DIRECTORY is empty.");
            System.exit(1);
        }
        theDir = new File(WORKING_DIRECTORY);
        if (theDir.exists()) {
            System.out.println("A pangenome database already exists at " + WORKING_DIRECTORY + GRAPH_DATABASE_PATH + ".");
            System.out.println("Do you want to connect to it? otherwisw it would be removed [y/n]? ");
            str = s.nextLine().toLowerCase();
            while (!str.equals("y") && !str.equals("n")){
                 System.out.println("Do you want to connect to it? otherwisw it would be removed [y/n]? ");
                 str = s.nextLine().toLowerCase();
            }
            if (str.equals("y")) {
                connect_panproteome();
                return;
            } else {
                try {
                    FileUtils.deleteRecursively(new File(WORKING_DIRECTORY + GRAPH_DATABASE_PATH));
                } catch (IOException ioe) {
                    System.out.println("Failed to delete the graph database");
                    System.exit(1);  
                }
            }
        }
        try {
            theDir.mkdir();
        } catch (SecurityException se) {
            System.out.println("Failed to create directory " + WORKING_DIRECTORY);
            System.exit(1);
        }
        graphDb = new GraphDatabaseFactory().newEmbeddedDatabaseBuilder(new File(WORKING_DIRECTORY + GRAPH_DATABASE_PATH))
            .setConfig(keep_logical_logs, "4 files").newGraphDatabase();  
        registerShutdownHook(graphDb);
    }

    /**
     * Connects to graph databases of the pan-genome.
     */
    public static void connect_panproteome(){
        if (new File(WORKING_DIRECTORY).exists()){
            if (! new File(WORKING_DIRECTORY + GRAPH_DATABASE_PATH).exists()) {
                System.out.println("No graph database found at " + WORKING_DIRECTORY);
                System.exit(1);
            }
            if (! new File(WORKING_DIRECTORY + INDEX_DATABASE_PATH).exists()) {
                System.out.println("No index database found at " + WORKING_DIRECTORY);
                System.exit(1);
            }
            graphDb = new GraphDatabaseFactory().newEmbeddedDatabaseBuilder(
                 new File(WORKING_DIRECTORY + GRAPH_DATABASE_PATH))
                .setConfig(keep_logical_logs, "4 files").newGraphDatabase();  
            registerShutdownHook(graphDb);
        } else {
            System.out.println("No panproteome found at " + WORKING_DIRECTORY);
            System.exit(1);
        }
    }    
    
    /**
     * Disconnects graph databases of the pan-genome.
     */
    public static void disconnect_panproteome(){
        graphDb.shutdown();
        File directory = new File(WORKING_DIRECTORY + GRAPH_DATABASE_PATH);
        for (File f : directory.listFiles()) {
            if (f.getName().startsWith("neostore.transaction.db.")) {
                f.delete();
            }
        }
    }    

    /**
     * Rebuilds the genome database from the graph database.
     */
    public static void rebuild_genome_database(){
    // read genomes information from the graph and rebuild the genomes database
        int genome, seqience, begin, end, j, len;
        long byte_number = 0;
        genomeDb = new SequenceDatabase(WORKING_DIRECTORY + GENOME_DATABASE_PATH, graphDb);
        StringBuilder seq = new StringBuilder();
        for (genome = 1; genome <= genomeDb.num_genomes; ++genome) {
            for (seqience = 1; seqience <= genomeDb.num_sequences[genome]; ++seqience) {
                begin = 1;
                end = (int) genomeDb.sequence_length[genome][seqience];
                extract_sequence_from_graph(seq, genome, seqience, begin, end);
                len = seq.length();
                if (len % 2 == 1) {
                    --len;
                }
                for (j = 0; j < len; j += 2, ++byte_number) {
                    genomeDb.genomes_buff[(int) (byte_number / genomeDb.MAX_BYTE_COUNT)].put((byte) ((binary[seq.charAt(j)] << 4) | binary[seq.charAt(j + 1)]));
                }
                if (len == seq.length() - 1) {
                    genomeDb.genomes_buff[(int) (byte_number / genomeDb.MAX_BYTE_COUNT)].put((byte) (binary[seq.charAt(len)] << 4));
                    ++byte_number;
                }
            }
        }
    }

    /**
     * Extracts a genomic region form the graph database.
     * 
     * @param seq The StringBuilder to write the sequence in. 
     * @param genome The query genome
     * @param sequence The query sequence
     * @param begin The start of the region
     * @param end The end of the region
     */
    public static void extract_sequence_from_graph(StringBuilder seq, int genome, int sequence, int begin, int end) {
        Relationship rel;
        Node neighbor, node;
        IndexPointer start_ptr;
        int loc, len = 0, node_len, neighbor_len, seq_len, position;
        String rel_name, origin;
        --begin;
        --end;
        origin = "G" + genome + "S" + sequence;
        seq_len = end - begin + 1;
        seq.setLength(0);
        start_ptr = locate(graphDb, genomeSc, indexSc, genome, sequence, begin);
        position = start_ptr.offset;
        node = graphDb.getNodeById(start_ptr.node_id);
        node_len = (int) node.getProperty("length");
    // Takes the part of the region lies in the first node of the path that region takes in the graph    
        if (start_ptr.canonical) {
            if (position + seq_len - 1 <= node_len - 1) { // The whole sequence lies in this node
                len += append_fwd(seq, (String) node.getProperty("sequence"), position, position + seq_len - 1);
            } else {
                len += append_fwd(seq, (String) node.getProperty("sequence"), position, node_len - 1);
            }
        } else {
            if (position - (seq_len - 1) >= 0) { // The whole sequence lies in this node
                len += append_rev(seq, (String) node.getProperty("sequence"), position - (seq_len - 1), position);
            } else {
                len += append_rev(seq, (String) node.getProperty("sequence"), 0, position);
            }
        }
    //  traverse the path of the region   
        while (len < seq_len) {
            //System.out.println(node.getId()+" "+len + " " + seq_len);
            loc = (begin + len) - K_SIZE + 1;
            rel = get_outgoing_edge(node, origin, loc);
            neighbor = rel.getEndNode();
            rel_name = rel.getType().name();
            neighbor_len = (int) neighbor.getProperty("length");
            if (rel_name.charAt(1) == 'F') {// Enterring forward side
                if (len + neighbor_len - K_SIZE + 1 > seq_len) // neighbor is the last node of the path
                    len += append_fwd(seq, (String) neighbor.getProperty("sequence"), K_SIZE - 1, seq_len - len + K_SIZE - 2);
                else 
                    len += append_fwd(seq, (String) neighbor.getProperty("sequence"), K_SIZE - 1, neighbor_len - 1);
            }else{ // Enterring reverse side
                if (len + neighbor_len - K_SIZE + 1 > seq_len) // neighbor is the last node of the pat
                    len += append_rev(seq, (String) neighbor.getProperty("sequence"), neighbor_len - K_SIZE - (seq_len - len) + 1, neighbor_len - K_SIZE);
                else 
                    len += append_rev(seq, (String) neighbor.getProperty("sequence"), 0, neighbor_len - K_SIZE);
            }
            node = neighbor;
        } // while
    }
    
    /**
     * Prints the manual of the software.
     */
    private static void print_help_message() {
        System.out.println("PanTools version 2.0\n" +
                "Welcome, Lygeri, to this new and improved version of Pantools! And may the codes be ever in your favor.\n" +
                "Get it? Codes? Because... odds... ha ha ha.\n" +
"PanTools is pan-genomic toolkit for comparative analysis of large number of genomes developed in Bioinformatics group of Wageningen University and Research Renter, the Netherlands. Please cite the relevant publication(s) from the list of publications if you use PanTools in your research.\n" +
"\n" +
"Licence\n" +
"PanTools has been licensed under GNU GENERAL PUBLIC LICENSE version 3.\n" +
"\n" +
"Publications\n" +
"\n" +
"PanTools: representation, storage and exploration of pan-genomic data.\n" +
"Efficient inference of homologs in large eukaryotic pan-proteomes.\n" +
"Pan-genomic read mapping\n" +
"\n" +
"\n" +
"Functionalities\n" +
"PanTools currently provides these functionalities:\n" +
"\n" +
"Construction of pan-genome\n" +
"Construction of pan-proteome\n" +
"Adding new genomes to the pan-genome\n" +
"Adding strucureal annotaions to the genomes\n" +
"Detecting homology groups based on similarity of proteins\n" +
"Retrieving features/regions/genomes\n" +
"Read mappping\n" +
"\n" +
"\n" +
"Requirements\n" +
"\n" +
"\n" +
"Java Virtual Machine version 1.8 or higher,\n" +
"Add path to the java executable to your OS path environment variable.\n" +
"\n" +
"KMC: A disk-based k-mer counter,\n" +
"After downloading the appropriate version (linux, macos or windows), add path to the kmc and kmc_tools executables to your OS path environment variable.\n" +
"\n" +
"MCL: The Markov Clustering Algorithm,\n" +
"After downloading and compiling the software, add path to the mcl executable to your OS path environment variable.\n" +
"\n" +
"\n" +
"Running the program\n" +
"Add path to the java archive of PanTools, located in the /dist sub-directory of PanTools project, to our OS path environment variable. Then run PanTools from command line by:\n" +
"$java <JVM options> -jar pantools.jar <sub-command> <arguments>\n" +
"\n" +
"Useful JVM options\n" +
"\n" +
"\n" +
"-server : To optimize JIT compilations for higher performance\n" +
"\n" +
"-Xmn(a number followed by m/g) : Minimum heap size in mega/giga bytes\n" +
"\n" +
"-Xmx(a number followed by m/g) : Maximum heap size in mega/giga bytes\n" +
"\n" +
"\n" +
"Sub-commands\n" +
"\n" +
"\n" +
"build_pangenome or bpg : To build a pan-genome out of a set of genomes.\n" +
"arguments:\n" +
"\n" +
"\n" +
"--database_path or -dp : Path to the pangenome database.\n" +
"\n" +
"--genomes-file or -gf : A text file containing paths to FASTA files of genomes;       each in a seperated line.\n" +
"\n" +
"--kmer-size or -ks : The size of k-mers. If it is not given or is out of valid range (6 <= K_SIZE <= 255), then an optimal value would be calculated automatically.\n" +
"\n" +
"\n" +
"\n" +
"build_panproteome or bpp : To build a pan-proteome out of a set of proteins.\n" +
"arguments:\n" +
"\n" +
"\n" +
"--database_path or -dp : Path to the pangenome database.\n" +
"\n" +
"--proteomes_file or -pf : A text file containing paths to FASTA files of proteomes; each in a seperated line.\n" +
"\n" +
"\n" +
"\n" +
"add_genomes or ag : To add new genomes to an available pan-genome.\n" +
"arguments:\n" +
"\n" +
"\n" +
"--database_path or -dp : Path to the pangenome database.\n" +
"\n" +
"--genomes-file or -gf : A text file containing paths to FASTA files of the new genomes to be added to the pangeome; each in a seperated line.\n" +
"\n" +
"\n" +
"\n" +
"add_annotations or aa : To add annotations to existing genomes.\n" +
"arguments:\n" +
"\n" +
"\n" +
"--database_path or -dp : Path to the pangenome database.\n" +
"\n" +
"--output-path or -op (default value: Database path determined by -dp) : Path to the output files.\n" +
"\n" +
"--annotations-file or -af : A text file each line of which contains genome number and path to the corresponding GFF file seperated by one space. Genomes are numbered in the same order they have been added to the pangenome. The protein sequence of the annotated genes will be also stored in the folder \"/proteins\" in the output path.\n" +
"\n" +
"--connect_annotations or -ca : Connects the annotated genomic features to the nodes of gDBG.\n" +
"\n" +
"\n" +
"\n" +
"retrieve_features or rf : To retrieve the sequence of annotated features from the pan-genome. For each genome a FASTA file containing the retrieved features will be stored in the output path. For example, genes.1.fasta contains all the genes annotated in genome 1.\n" +
"arguments:\n" +
"\n" +
"\n" +
"--database_path or -dp : Path to the pangenome database.\n" +
"\n" +
"--output-path or -op (default value: Database path determined by -dp) : Path to the output files.\n" +
"\n" +
"--genome-numbers or -gn : A text file containing genome_numbers for which the features will be retrieved.\n" +
"\n" +
"\n" +
"\n" +
"--feature-type or -ft (default value: gene) : The feature name; for example gene, mRNA, exon, tRNA, etc.\n" +
"\n" +
"\n" +
"retrieve_regions or rr : To retrieve the sequence of some genomic regios from the pan-genome. The resulting FASTA files will be stored in the output path.\n" +
"arguments:\n" +
"\n" +
"\n" +
"--database_path or -dp : Path to the pangenome database.\n" +
"\n" +
"--output-path or -op (default value: Database path determined by -dp) : Path to the output files.\n" +
"\n" +
"--regions-file or -rf : A text file containing records with genome_number, sequence_number, begin and end positions seperated by one space for each region.\n" +
"\n" +
"\n" +
"\n" +
"retrieve_genomes or rg : To retrieve the full sequence of some genomes. The resulting FASTA files will be stored in the output path.\n" +
"arguments:\n" +
"\n" +
"\n" +
"--database_path or -dp : Path to the pangenome database.\n" +
"\n" +
"--output-path or -op (default value: Database path determined by -dp) : Path to the output files.\n" +
"\n" +
"--genome-numbers or -gn : A text file containing genome_numbers to be retrieved in each line..\n" +
"\n" +
"\n" +
"\n" +
"group or g : To create homology groups in the protein space of the pan-genome (pan-proteome). The resulting homology groups will be stored in the output path.\n" +
"arguments:\n" +
"\n" +
"\n" +
"--database_path or -dp : Path to the pangenome database.\n" +
"\n" +
"--output-path or -op (default value: Database path determined by -dp) : Path to the output files.\n" +
"\n" +
"--intersection-rate or -ir (default valuue: 0.09, valid range: [0.001..0.1]) : The fraction of k-mers needs to be shared by two intersecting proteins.\n" +
"\n" +
"--min-normalized-similarity or -mns (default value: 95, valid range: [1..99]) : The minimum normalized similarity score of two proteins.\n" +
"\n" +
"--mcl-inflation or -mi (default value: 9.6, valid range: (1..19)): The MCL inflation.\n" +
"\n" +
"--contrast or -ct (default value: 8, valid range: (0..10)) : The contrast factor.\n" +
"\n" +
"--relaxation or rn (default value: 1, valid range: [1..8]) : The relaxation in homology calls.\n" +
"\n" +
"--threads-number or -tn (default value: 1) : The number of parallel working threads.\n" +
"\n" +
"\n" +
"\n" +
"map or m : To map single or paired-end short reads to all or a sebset of the constituent genomes. The resulting SAM/BAM files will be stored in the output path.\n" +
"argument:\n" +
"\n" +
"\n" +
"--database_path or -dp : Path to the pangenome database.\n" +
"\n" +
"\n" +
"-1 : The first short-read archive in FASTQ format, which can be gz/bz2 compressed. This file can be precessed interleaved by -il option.\n" +
"\n" +
"\n" +
"-2 : Optionally, the second short-read archive in FASTQ format, which can be gz/bz2 compressed.\n" +
"\n" +
"\n" +
"--genome-numbers or -gn : A text file containing genome_numbers to map reads against in each line.\n" +
"\n" +
"\n" +
"--output-path or -op (default value: Database path determined by -dp) : Path to the output files.\n" +
"\n" +
"\n" +
"--threads-number or -tn (default value: 1) : The number of parallel working threads.\n" +
"\n" +
"\n" +
"--min-mapping-identity or -mmi (default value: 0.5, valid range: [0..1)) : The minimum acceptable identity of the alignment.\n" +
"\n" +
"\n" +
"--num-kmer-samples or *-nks (default value: 15, valid range: [1..r-k+1]) : The number of kmers sampled from read.\n" +
"\n" +
"\n" +
"--min-hit-length or -mhl (default value: 13, valid range: [10..100]) : The minimum acceptable length of alignment after soft-clipping.\n" +
"\n" +
"\n" +
"--max-alignment-length or -mal (default value: 1000, valid range: [50..5000]) : The maximum acceptable length of alignment.\n" +
"\n" +
"\n" +
"--max-fragment-length or -mfl (default value: 2000, valid range: [50..5000]) : The maximum acceptable length of fragment.\n" +
"\n" +
"\n" +
"--max-num-locations or -mnl (default value: 15, valid range: [1..100]) : The maximum number of location of candidate hits to examine.\n" +
"\n" +
"\n" +
"--alignment-band or -ab (default value: 5, valid range: [1..100]) : The length of bound of banded alignment.\n" +
"\n" +
"\n" +
"--clipping-stringency or -ci (default value: 1) : The stringency of soft-clipping.\n" +
"0 : no soft clipping\n" +
"1 : low\n" +
"2 : medium\n" +
"3 : high\n" +
"\n" +
"\n" +
"--bam-format or -bf : Writes the alignment files in .BAM format.\n" +
"\n" +
"\n" +
"--alignment-mode or -am (default value: 2) : The alignment mode.\n" +
"-1 : Competitive, none-bests\n" +
"-2 : Competitive, random-best\n" +
"-3 : Competitive, all-bests\n" +
"1 : Normal, none-bests\n" +
"2 : Normal, random-best\n" +
"3 : Normal, all-bests\n" +
"0 : Normal, all-hits\n" +
"\n" +
"\n" +
"--interleaved or -i : Process the fastq file as an interleaved paired-end archive.\n" +
"\n" +
"\n" +
"\n" +
"\n" +
"version or v : To show the versions of PanTools, JVM and Neo4j.\n" +
"\n" +
"\n" +
"help or h: To show the mannual of the tool.\n" +
"\n" +
"\n" +
"\n" +
"Visualization in the Neo4j browser\n" +
"Neo4j browser allows you to visualize parts of the pan-genome graph and run Cypher queries and receive the results in a tabular or a graphic format. You need to download the appropriate version of Neo4j (use version sub-command to see the consistent version). To visualize a pan-genome:\n" +
"\n" +
"Add path to the Neo4j /bin directory to the path environment variable.\n" +
"Add path to your pan-genome in the Neo4j configuration file NEO4J-DIRECTORY/conf/neo4j.conf:\n" +
"dbms.directories.data = PATH_TO_THE_PANGENOME_DATABASE\n" +
"Start the Neo4j database server from the shell:\n" +
"\n" +
"   $neo4j start\n" +
"\n" +
"open an internet browser and open http://localhost:7474.\n" +
"To visualize the whole pangenome, type this simple Cypher command in the browser:\n" +
"MATCH (n) RETURN n\n" +
"Stop the Neo4j database server from the shell:\n" +
"\n" +
"$neo4j stop");
    }

    /**
     * Estimates and prints the peak memory used during the execution of the program. 
     */
    public static void print_peak_memory() {
        long memoryUsage = 0;
        try {
            for (MemoryPoolMXBean pool : ManagementFactory.getMemoryPoolMXBeans()) {
                MemoryUsage peak = pool.getPeakUsage();
                memoryUsage += peak.getUsed();
            }
            System.out.println("Peak memory : " + memoryUsage / 1024 / 1024 + " MB");
        } catch (Throwable t) {
            System.err.println("Exception in agent: " + t);
        }
    }
    
    /**
     * Writes a sequence in a FASTA file with specified length for lines.
     * 
     * @param fasta_file The FASTA file object
     * @param seq The sequence to be written in the file.
     * @param length Length of the lines.
     */    
    public static void write_fasta(BufferedWriter fasta_file, String seq, int length) {
        int i;
        try {
            for (i = 1; i <= seq.length(); ++i) {
                fasta_file.write(seq.charAt(i - 1));
                if (i % length == 0) {
                    fasta_file.write("\n");
                }
            }
            fasta_file.write("\n");
        } catch (IOException ioe) {

        }

    }    
    
    /**
     * Calsulates the reverse complement of a given string.
     * 
     * @param s The input string
     * @return The reverse complement of the input string
     */     
    public static void reverse_complement(StringBuilder s) {
        char ch;
        int i, j;
        for (i = 0, j = s.length() - 1; i < j; ++i, --j) {
            ch = s.charAt(i);
            s.setCharAt(i, complement(s.charAt(j)));
            s.setCharAt(j, complement(ch));
        }
        if (i == j)
            s.setCharAt(i, complement(s.charAt(i)));
    }

    /**
     * Calculates the complement of an IUPAC symbol
     * 
     * @param ch The input symbol
     * @return The complement symbol
     */
    public static char complement(char ch) {
        switch (ch) {
            case 'A':
                return 'T';
            case 'C':
                return 'G';
            case 'G':
                return 'C';
            case 'T':
                return 'A';
            case 'R':
                return 'Y';
            case 'Y':
                return 'R';
            case 'K':
                return 'M';
            case 'M':
                return 'K';
            case 'B':
                return 'V';
            case 'V':
                return 'B';
            case 'D':
                return 'H';
            case 'H':
                return 'D';
            default:
                return ch;
        }
    }
    
    /**
     * Executes a shell command. 
     * 
     * @param command The command
     * @return The output of the bash command
     */
    public static String executeCommand(String command) {
        StringBuilder exe_output = new StringBuilder();
        String line = "";
        Process p;
        try {
            p = Runtime.getRuntime().exec(command);
            p.waitFor();
            BufferedReader reader = new BufferedReader(new InputStreamReader(p.getInputStream()));
            while ((line = reader.readLine()) != null) {
                exe_output.append(line + "\n");
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        return exe_output.toString();
    }    

    /**
     * Executes a shell command in a limited number of seconds.
     * 
     * @param command The command
     * @param seconds The number of seconds
     * @return The output of the bash command
     */
    public static boolean executeCommand_for(String command, int seconds) {
        Process p;
        boolean success = false;
        try {
            p = Runtime.getRuntime().exec(command);
            success = p.waitFor(seconds, TimeUnit.SECONDS);
        } catch (Exception e) {
            e.printStackTrace();
        }
        return success;
    }
    
    /**
     * Determines if a file is in FASTA format.
     * 
     * @param file_name Path to the file
     * @return True if file is in FASTA format, or False
     */
    public static boolean is_fasta(String file_name){
        try{
            BufferedReader in = open_file(file_name);
            String line;
            while ((line = in.readLine()) != null){
                line = line.trim();
                if (line.equals("")) 
                    continue;
                else {
                    in.close();
                    return line.charAt(0) == '>';
                }            
            }
        } catch (IOException ex){
            System.err.println(ex.getMessage());
        }
        return false;
    }

    /**
     * Determines if a file is in FASTQ format.
     * 
     * @param file_name Path to the file
     * @return True if file is in FASTQ format, or False
     */
    public static boolean is_fastq(String file_name){
        try{
            BufferedReader in = open_file(file_name);
            String line;
            while ((line = in.readLine()) != null){
                line = line.trim();
                if (line.equals("")) 
                    continue;
                else {
                    in.close();
                    return line.charAt(0) == '@';
                }
            }
        } catch (IOException ex){
            System.err.println(ex.getMessage());
        }
        return false;
    }
    
    /**
     * Opens a possibly compressed file.
     * 
     * @param filename Path to the file
     * @return The buffered reader to the input file 
     */
    public static BufferedReader open_file(String filename){
        try{        
            String[] fields = filename.split("\\.");
            String file_type = fields[fields.length - 1].toLowerCase();
            if (file_type.equals("gz") || file_type.equals("gzip") || file_type.equals("bz2") || file_type.equals("bzip2"))
                return getBufferedReaderForCompressedFile(filename);//BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(filename)), "UTF-8"));                    
            else 
                return new BufferedReader(new BufferedReader(new FileReader(filename)));                    
        } catch (IOException ex){
            System.out.println(ex.getMessage());
            return null;
        }
    }
    
    /**
     * Counts the number of lines in a file.
     * 
     * @param file_name Path to the file
     * @param skip_empty determines if empty lines should be skipped
     * @return The count of the lines
     */
    public static int get_lines_count(String file_name, boolean skip_empty){
        int count = 0;
        try{
            BufferedReader in = open_file(file_name);
            String line;
            while ((line = in.readLine()) != null){
                line = line.trim();
                if (skip_empty && line.equals("")) 
                    continue;
                ++count;
            }
            in.close();
        } catch (IOException ex){
            System.err.println(ex.getMessage());
        }
        return count;
    }
    
    /**
     * Gives the buffered reader object for a compressed file.
     * 
     * @param fileIn Path of the input file
     * @return 
     */
    public static BufferedReader getBufferedReaderForCompressedFile(String fileIn){
        try{
            FileInputStream fin = new FileInputStream(fileIn);
            BufferedInputStream bis = new BufferedInputStream(fin);
            CompressorInputStream input = new CompressorStreamFactory().createCompressorInputStream(bis);
            BufferedReader br2 = new BufferedReader(new InputStreamReader(input));
            return br2;
        } catch (Exception ex){
            System.err.println(ex.getMessage() + "\nFailed to open the compresse file!");
            return null;
        }
    }    

    /**
     * Shuts down the graph database if the program halts unexpectedly.
     * 
     * @param graphDb The graph database object 
     */
    public static void registerShutdownHook(final GraphDatabaseService graphDb) {
        Runtime.getRuntime().addShutdownHook(new Thread() {
            @Override
            public void run() {
                graphDb.shutdown();
            }
        });
    }

     /*
    #NEW
    */
    public static void check_database(){ // NEW
        File file = new File(PATH_TO_THE_PANGENOME_DATABASE);
        if (!file.exists()) {
            System.out.println("\nThe provided database was not found!\n");
            System.exit(1);
        } else if (!new File(PATH_TO_THE_PANGENOME_DATABASE + GRAPH_DATABASE_PATH).exists()) {
            System.out.println("\nUnable to open the database provided via -dp");
            System.exit(0);
        }
    }

    private static void print_eeftools(){
        System.out.println("Welcome back... Here we go again!\n" +
            " ▄▄▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄            ▄▄▄▄▄▄▄▄▄▄▄ \n" +
            "▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░▌          ▐░░░░░░░░░░░▌\n" +
            "▐░█▀▀▀▀▀▀▀▀▀ ▐░█▀▀▀▀▀▀▀▀▀ ▐░█▀▀▀▀▀▀▀▀▀  ▀▀▀▀█░█▀▀▀▀ ▐░█▀▀▀▀▀▀▀█░▌▐░█▀▀▀▀▀▀▀█░▌▐░▌          ▐░█▀▀▀▀▀▀▀▀▀ \n" +
            "▐░▌          ▐░▌          ▐░▌               ▐░▌     ▐░▌       ▐░▌▐░▌       ▐░▌▐░▌          ▐░▌          \n" +
            "▐░█▄▄▄▄▄▄▄▄▄ ▐░█▄▄▄▄▄▄▄▄▄ ▐░█▄▄▄▄▄▄▄▄▄      ▐░▌     ▐░▌       ▐░▌▐░▌       ▐░▌▐░▌          ▐░█▄▄▄▄▄▄▄▄▄ \n" +
            "▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌     ▐░▌     ▐░▌       ▐░▌▐░▌       ▐░▌▐░▌          ▐░░░░░░░░░░░▌\n" +
            "▐░█▀▀▀▀▀▀▀▀▀ ▐░█▀▀▀▀▀▀▀▀▀ ▐░█▀▀▀▀▀▀▀▀▀      ▐░▌     ▐░▌       ▐░▌▐░▌       ▐░▌▐░▌           ▀▀▀▀▀▀▀▀▀█░▌\n" +
            "▐░▌          ▐░▌          ▐░▌               ▐░▌     ▐░▌       ▐░▌▐░▌       ▐░▌▐░▌                    ▐░▌\n" +
            "▐░█▄▄▄▄▄▄▄▄▄ ▐░█▄▄▄▄▄▄▄▄▄ ▐░▌               ▐░▌     ▐░█▄▄▄▄▄▄▄█░▌▐░█▄▄▄▄▄▄▄█░▌▐░█▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄█░▌\n" +
            "▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░▌               ▐░▌     ▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌\n " +
            "▀▀▀▀▀▀▀▀▀▀▀  ▀▀▀▀▀▀▀▀▀▀▀  ▀                 ▀       ▀▀▀▀▀▀▀▀▀▀▀  ▀▀▀▀▀▀▀▀▀▀▀  ▀▀▀▀▀▀▀▀▀▀▀  ▀▀▀▀▀▀▀▀▀▀▀\n"
                + "--------------------------------------------          -------------------------------------------");
    }
       
   
    
    
    
    static class draad extends Thread{
        private String naam;
        private int aantal;
        private int kolom;
        private String s;
  
        public draad(String threadName, String naam, int aantal, int kolom) { // constructor
            super(threadName);
            this.naam = naam;
            this.aantal = aantal;
            this.kolom = kolom;
     
            s = " ";
            for (int j = 0; j < kolom; j++)
                s = "\t" + s;
                s = s + "x";
        }
  
        public void run() {
            for (int i = 0; i < aantal; i++) {
                System.out.println(naam + i + s);
            }
                //System.out.println("hoi" + kolom);
            EefWP1.write_string_to_file_here("", "test" + kolom);
        }         
    }
     
    public static void create_full_path_working_directory(){
        //System.out.println(current_path + "\n" + WORKING_DIRECTORY);
        int count = StringUtils.countMatches(WORKING_DIRECTORY, "../"); 
        String [] path_array = current_path.split("/");
        String temp = WORKING_DIRECTORY.replace ("../","");
        WD_full_path = "";
        for (int j=0; j< path_array.length-count; j++) {
            WD_full_path += path_array[j] + "/";
        }
        WD_full_path  += temp;
        //System.out.println("dit: " + WD_full_path );
        if (!new File(WD_full_path).exists()){      
            WD_full_path =  WORKING_DIRECTORY;
        }
    }
}
