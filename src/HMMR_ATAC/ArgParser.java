package HMMR_ATAC;

import java.io.File;

public class ArgParser {

    private String[] args;

    //Required inputs
    private File bam = null;
    private File index = null;
    private String genomeFile = null;
    //private  String bigWig = null;

    //Optional Inputs
    private String means;//comma separted list of intial mean values for frag dist
    private String stddevs;//comma separted list of initial standard deviations for frag dist
    private boolean fragEM = true; //whether or not to perform fragment dist em
    private int minMapQ = 30; //minimum mapping quality of reads to keep
    private int lower = 10; //lower bound for fold change range for choosing training sites
    private int upper = 20; //upper bound for fold change range for choosing training sites
    private float zscore = 100; //zscored read coverage to exclude from viterbi decoding
    private float zscoreprescan = 0; //above this zscore to be included in viterbi decoding
    private String output; //output name 
    private boolean peaks = true; // whether to print peaks
    private boolean bg = false; // whether to print bedgraph
    private String blacklist = null;
    private boolean keepDups = false;
    private int minLength = 200;
    private String scoreSys = "max";
    private boolean bgScore = false;
    private int k = 3;
    private int trim = 0;
    private String trainingRegions;
    private int vitWindow = 25000000;
    private File modelFile = null;
    private boolean stopAfterModel = false;

    /**
     * Main constructor
     *
     * @param an Array of Strings passed from the main program
     */
    public ArgParser(String[] a) {
        args = a;
        set();
    }

    /**
     * Access stop after model
     *
     * @return a boolean to determine whether to stop the program after model
     * generation
     */
    public boolean getStopAfterModel() {
        return stopAfterModel;
    }

    /**
     * Access score the bedgraph
     *
     * @return a boolean to determine whether to score the bedgraph
     */
    public boolean getBGScore() {
        return bgScore;
    }

    /**
     * Access the score system
     *
     * @return a String representing the peak score system
     */
    public String getScore() {
        return scoreSys;
    }

    /**
     * Access the input BAM file
     *
     * @return a File representing the BAM input file
     */
    public File getBam() {
        return bam;
    }

    /**
     * Access the BAM index file
     *
     * @return a File representing the BAM index file
     */
    public File getIndex() {
        return index;
    }

    /**
     * Access the genome file
     *
     * @return a String representing the genome file
     */
    public String getGenome() {
        return genomeFile;
    }

    /**
     * Access the bigwig file
     *
     * @return a String representing the bigwig file
     */
//		public String getBigWig(){return bigWig;}
    /**
     * Access the means of the mixture model
     *
     * @return a String representing the means of the mixture model
     */
    public String getMeans() {
        return means;
    }

    /**
     * Access the means of the mixture model
     *
     * @return a String representing the standard deviation of the mixture model
     */
    public String getStd() {
        return stddevs;
    }

    /**
     * Access EM training
     *
     * @return a boolean to determine whether to perform EM training
     */
    public boolean getEM() {
        return fragEM;
    }

    /**
     * Access minimum mapping quality
     *
     * @return an integer representing the minimum mapping quality score of the
     * reads to use in track generation
     */
    public int getMinQ() {
        return minMapQ;
    }

    /**
     * Access the lower bounds of fold change to identify training regions
     *
     * @return an interger representing the lower bounds
     */
    public int getLower() {
        return lower;
    }

    /**
     * Access the upper bounds of fold change to identify training regions
     *
     * @return an interger representing the upper bounds
     */
    public int getUpper() {
        return upper;
    }

    /**
     * Access the zscored value of regions to exclude from training or decoding
     *
     * @return an interger representing the zscore maximum
     */
    public float getZscore() {
        return zscore;
    }

    /**
     * Access the minimum zscored value of regions to include from training or decoding
     *
     * @return an interger representing the zscore minimum
     */
    public float getZscorePrescan() {
        return zscoreprescan;
    }
    
    /**
     * Access the output prefix name
     *
     * @return a String representing the output file names
     */
    public String getOutput() {
        return output;
    }

    /**
     * Access report peaks
     *
     * @return a boolean to determine whether to report peaks
     */
    public boolean getPeaks() {
        return peaks;
    }

    /**
     * Access report bedgraph
     *
     * @return a boolean to determine whether to report genome wide bedgraph
     */
    public boolean getBedgraph() {
        return bg;
    }

    /**
     * Access the BED file name for excluded regions
     *
     * @return a String representing the name of the BED file of regions to
     * exclude
     */
    public String getBlacklist() {
        return blacklist;
    }

    /**
     * Access keep duplicates
     *
     * @return a boolean to determine whether to keep duplicate reads
     */
    public boolean getKeepDups() {
        return keepDups;
    }

    /**
     * Access minimum length for peaks to report
     *
     * @return an integer representing the minimum length of peaks to be
     * reported
     */
    public int getMinLength() {
        return minLength;
    }

    /**
     * Access the number of states in the model
     *
     * @return an integer representing the number of states in the model
     */
    public int getK() {
        return k;
    }

    public int getTrim() {
        return trim;
    }

    /**
     * Access the number of distributions to trim
     *
     * @return an integer representing the number of distributions to trim off
     * the data tracks
     */
    public String getTrainingRegions() {
        return trainingRegions;
    }

    /**
     * Access the size of the region to decode with viterbi algorithm
     *
     * @return an integer representing the size, in bp, of the genomic regions
     * to run viterbi on
     */
    public int getWindow() {
        return vitWindow;
    }

    /**
     * Access the inputed model file
     *
     * @return a File representing the binary model file generated from a
     * previous HMMR run
     */
    public File getModelFile() {
        return modelFile;
    }

    /**
     * Parse the argument string and set the variables
     */
    private void set() {
        for (int i = 0; i < args.length; i++) {

            switch (args[i].charAt(0)) {
                case '-':

                    switch (args[i].charAt((1))) {

                        case 'b':
                            bam = new File(args[i + 1]);
                            i++;
                            break;
                        case 'i':
                            index = new File(args[i + 1]);
                            i++;
                            break;
                        case 'g':
                            genomeFile = args[i + 1];
                            i++;
                            break;
//				case'w':
//					bigWig = args[i+1];
//					i++;
//					break;
                        case 'm':
                            means = args[i + 1];
                            i++;
                            break;
                        case 's':
                            stddevs = args[i + 1];
                            i++;
                            break;
                        case 'f':
                            String temp = args[i + 1].toLowerCase();
                            if (temp.contains("t")) {
                                fragEM = true;
                            } else {
                                fragEM = false;
                            }
                            i++;
                            break;
                        case 'q':
                            minMapQ = Integer.parseInt(args[i + 1]);
                            i++;
                            break;
                        case 'u':
                            upper = Integer.parseInt(args[i + 1]);
                            i++;
                            break;
                        case 'l':
                            lower = Integer.parseInt(args[i + 1]);
                            i++;
                            break;
                        case 'z':
                            zscore = Float.parseFloat(args[i + 1]);
                            i++;
                            break;
                         case 'Z':
                            zscoreprescan = Float.parseFloat(args[i + 1]);
                            i++;
                            break;
                        case 'o':
                            output = args[i + 1];
                            i++;
                            break;
                        case 'e':
                            blacklist = args[i + 1];
                            i++;
                            break;
                        case 'p':
                            String t = args[i + 1].toLowerCase();
                            if (t.contains("t")) {
                                peaks = true;
                            } else {
                                peaks = false;
                            }
                            i++;
                            break;
                        case 'k':
                            k = Integer.parseInt(args[i + 1]);
                            i++;
                            break;
                        case 't':
                            trainingRegions = args[i + 1];
                            i++;
                            break;

                        case 'h':
                            printUsage();
                        //System.exit(0);
                        case '-':
                            switch (args[i].substring(2)) {
                                case "bam":
                                    bam = new File(args[i + 1]);
                                    i++;
                                    break;
                                case "index":
                                    index = new File(args[i + 1]);
                                    i++;
                                    break;
                                case "genome":
                                    genomeFile = args[i + 1];
                                    i++;
                                    break;
//					case"wig":
//						bigWig = args[i+1];
//						i++;
//						break;
                                case "means":
                                    means = args[i + 1];
                                    i++;
                                    break;
                                case "stddev":
                                    stddevs = args[i + 1];
                                    i++;
                                    break;
                                case "fragem":
                                    String t1 = args[i + 1].toLowerCase();
                                    if (t1.contains("t")) {
                                        fragEM = true;
                                    } else {
                                        fragEM = false;
                                    }
                                    i++;
                                    break;
                                case "minmapq":
                                    minMapQ = Integer.parseInt(args[i + 1]);
                                    i++;
                                    break;
                                case "upper":
                                    upper = Integer.parseInt(args[i + 1]);
                                    i++;
                                    break;
                                case "lower":
                                    lower = Integer.parseInt(args[i + 1]);
                                    i++;
                                    break;
                                case "zscore":
                                    zscore = Float.parseFloat(args[i + 1]);
                                    i++;
                                    break;
                                case "zscoreprescan":
                                    zscoreprescan = Float.parseFloat(args[i + 1]);
                                    i++;
                                    break;                                    
                                case "output":
                                    output = args[i + 1];
                                    i++;
                                    break;
                                case "blacklist":
                                    blacklist = args[i + 1];
                                    i++;
                                    break;
                                case "peaks":
                                    String t2 = args[i + 1].toLowerCase();
                                    if (t2.contains("t")) {
                                        peaks = true;
                                    } else {
                                        peaks = false;
                                    }
                                    i++;
                                    break;
                                case "bedgraph":
                                    String t3 = args[i + 1].toLowerCase();
                                    if (t3.contains("t")) {
                                        bg = true;
                                    } else {
                                        bg = false;
                                    }
                                    i++;
                                    break;
                                case "minlen":
                                    minLength = Integer.parseInt(args[i + 1]);
                                    i++;
                                    break;
                                case "score":
                                    scoreSys = args[i + 1].toLowerCase();
                                    i++;
                                    break;
                                case "bgscore":
                                    String t4 = args[i + 1].toLowerCase();
                                    if (t4.contains("t")) {
                                        bgScore = true;
                                    } else {
                                        bgScore = false;
                                    }
                                    i++;
                                    break;
                                case "kmeans":
                                    k = Integer.parseInt(args[i + 1]);
                                    i++;
                                    break;
                                case "trim":
                                    trim = Integer.parseInt(args[i + 1]);
                                    i++;
                                    break;
                                case "training":
                                    trainingRegions = args[i + 1];
                                    i++;
                                    break;
                                case "window":
                                    vitWindow = Integer.parseInt(args[i + 1]);
                                    i++;
                                    break;
                                case "model":
                                    modelFile = new File(args[i + 1]);
                                    i++;
                                    break;
                                case "modelonly":
                                    String TEMP = args[i + 1].toLowerCase();
                                    if (TEMP.contains("t")) {
                                        stopAfterModel = true;
                                    } else {
                                        stopAfterModel = false;
                                    }
                                    i++;
                                    break;
                                case "help":
                                    printUsage();
                                //System.exit(0);
                            }
                    }
            }
        }//for loop
    }

    /**
     * Print usage statement
     */
    public void printUsage() {
        System.out.println("Usage: java -jar HMMRATAC.jar");
        System.out.println("\nRequired Parameters:");
        System.out.println("\t-b , --bam <BAM> Sorted BAM file containing the ATAC-seq reads");
        System.out.println("\t-i , --index <BAI> Index file for the sorted BAM File");
        System.out.println("\t-g , --genome <GenomeFile> Two column, tab delimited file containing genome size stats");
//			System.out.println("\t-w , --wig <BigWig> Whole genome big wig file created using all reads");
        System.out.println("\nOptional Parameters:");
        System.out.println("\t-m , --means <double> Comma separated list of initial mean values for the fragment distribution. Default = 50,200,400,600");
        System.out.println("\t-s , --stddev <double> Comma separated list of initial standard deviation values for fragment distribution. Default = 20,20,20,20");
        System.out.println("\t-f , --fragem <true || false> Whether to perform EM training on the fragment distribution. Default = True");
        System.out.println("\t-q , --minmapq <int> Minimum mapping quality of reads to keep. Default = 30");
        System.out.println("\t-u , --upper <int> Upper limit on fold change range for choosing training sites. Default = 20");
        System.out.println("\t-l , --lower <int> Lower limit on fold change range for choosing training sites. Default = 10");
        System.out.println("\t-z , --zscore <float> Zscored read depth to mask during Viterbi decoding. Default = 100");
        System.out.println("\t-Z , --zscoreprescan <float> Minimum zscored read depth to be included in Viterbi decoding. Default = 0");
        System.out.println("\t-o , --output <File> Name for output files. Default = NA");
        System.out.println("\t-e , --blacklist <BED_File> bed file of blacklisted regions to exclude");
        System.out.println("\t-p , --peaks <true || false> Whether to report peaks in bed format. Default = true");
        System.out.println("\t-k , --kmeans <int> Number of States in the model. Default = 3. If not k=3, recommend NOT calling peaks, use bedgraph");
        System.out.println("\t-t , --training <BED_File> BED file of training regions to use for training model, instead of foldchange settings");
        System.out.println("\t--bedgraph <true || false> Whether to report whole genome bedgraph of all state anntations. Default = false");
        System.out.println("\t--minlen <int> Minimum length of open region to call peak. Note: -p , --peaks must be set. Default = 200");
        System.out.println("\t--score <max || ave || med || fc || zscore || all> What type of score system to use for peaks. Can be used for ranking peaks. Default = max");
        System.out.println("\t--bgscore <true || false> Whether to add the HMMR score to each state annotation in bedgraph. Note: this adds considerable time. Default = False");
        System.out.println("\t--trim <int> How many signals from the end to trim off (ie starting with tri signal then di etc). This may be useful if your data doesn't contain many large fragments. Default = 0");
        System.out.println("\t--window <int> Size of the bins to split the genome into for Viterbi decoding.\n\t To save memory, the genome is split into <int> long bins and viterbi decoding occurs across each bin. \n\tDefault = 25000000. Note: For machines with limited memory, it is recomended to reduce the size of the bins.");
        System.out.println("\t--model <File> Binary model file (generated from previous HMMR run) to use instead of creating new one");
        System.out.println("\t--modelonly <true || false> Whether or not to stop the program after generating model. Default = false");
        System.out.println("\t-h , --help Print this help message and exit.");
        System.exit(0);
    }

}
