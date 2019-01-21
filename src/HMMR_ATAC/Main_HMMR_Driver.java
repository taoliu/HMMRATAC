package HMMR_ATAC;

/*
 * Written by: Evan Tarbell 
 * evantarb@buffalo.edu
 */
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.logging.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.ObservationVector;
import be.ac.ulg.montefiore.run.jahmm.OpdfMultiGaussian;
import be.ac.ulg.montefiore.run.jahmm.io.HmmBinaryReader;
import be.ac.ulg.montefiore.run.jahmm.io.HmmBinaryWriter;
import ATACFragments.FragPileupGen;
import FormatConverters.PileupToBedGraph;
import GEMM.HMMR_EM;
import GenomeFileReaders.GenomeFileReader;
import GenomeFileReaders.bedFileReader;
import Node.PileupNode2;
import Node.ScoreNode;
import Node.TagNode;
import RobustHMM.KMeansToHMM;
import RobustHMM.RobustHMM;

//import WigMath.FoldChange;//10_13_18
//import WigMath.PullWigAboveCutoff;//10_13_18
import WigMath.bedGraphMath;
import WigMath.pileup;

public class Main_HMMR_Driver {

    //Required inputs
    private static File bam = null;
    private static File index = null;
    private static String genomeFile = null;
    //private static String bigWig = null;//10_13_18

    //Optional Inputs
    private static String means;//comma separted list of intial mean values for frag dist
    private static String stddevs;//comma separted list of initial standard deviations for frag dist
    private static boolean fragEM = true; //whether or not to perform fragment dist em
    private static int minMapQ = 30; //minimum mapping quality of reads to keep

    private static int lower = 10; //lower bound for fold change range for choosing training sites
    private static int upper = 20; //upper bound for fold change range for choosing training sites
    private static int zscore = 100; //zscored read coverage to exclude from viterbi decoding
    private static String output; //output name 
    private static boolean peaks; // whether to print peaks
    private static boolean bg; // whether to print bedgraph
    private static String blacklist = null;
    private static int minLength;
    private static String scoreSys;
    private static boolean BGScore;
    private static int k = 3;
    private static int trim = 0;
    private static int vitWindow = 25000000;
    private static File modelFile;
    private static boolean stopAfterModel = false;

    private static String trainingRegions;
    
    @SuppressWarnings("unchecked")
    public static void main(String[] args) throws IOException {

        ArgParser p = new ArgParser(args);
        bam = p.getBam();
        index = p.getIndex();
        genomeFile = p.getGenome();
        //bigWig=p.getBigWig();//10_13_18
        means = p.getMeans();
        stddevs = p.getStd();
        fragEM = p.getEM();
        minMapQ = p.getMinQ();
        lower = p.getLower();
        upper = p.getUpper();
        zscore = p.getZscore();
        output = p.getOutput();
        peaks = p.getPeaks();
        bg = p.getBedgraph();
        blacklist = p.getBlacklist();
        minLength = p.getMinLength();
        scoreSys = p.getScore();
        BGScore = p.getBGScore();
        k = p.getK();
        trim = p.getTrim();
        trainingRegions = p.getTrainingRegions();
        vitWindow = p.getWindow();
        modelFile = p.getModelFile();
        //For run time calculation
        Long startTime = System.currentTimeMillis();
        Long tmpTime = System.currentTimeMillis();

        //Declare output name
        if (output == null) {
            output = "NA";
        }
        
        //prepare logger
        Logger logger = Logger.getLogger("HMMRATAC") ;
        System.setProperty("java.util.logging.SimpleFormatter.format", "[%1$tc] %4$s: %5$s%n");

        //PrintStream log = new PrintStream(output + ".log");
        logger.info("HMMRATAC start...");

        //Exit program if BAM file or Index file not given
        if (bam == null || index == null || genomeFile == null) {//|| bigWig == null){//10_13_18
            p.printUsage();
            System.exit(1);
        }
        logger.config("Arguments Used:");
        for (int i = 0; i < args.length - 1; i += 2) {
            logger.config(args[i] + "\t" + args[i + 1]);
        }

        //Read in genome size stats
        GenomeFileReader gReader = new GenomeFileReader(genomeFile);
        ArrayList<TagNode> genomeStats = gReader.getMap();
        gReader = null;

        //Read in blacklisted if inputted
        ArrayList<TagNode> black = null;
        if (blacklist != null) {
            black = new bedFileReader(blacklist).getData();

        }

        /*
		 * Set fragment length distribution parameters. 
		 * Use inputed values to set initial values, if provided. 
		 * Else use defaults
         */
        double[] fragMeans = new double[4];
        double[] fragStddevs = new double[4];
        double[] mode = new double[4];
        mode[1] = mode[2] = mode[3] = 2;
        mode[0] = 0.5;

        if (means != null) {
            String[] mu = means.split(",");
            for (int i = 0; i < mu.length; i++) {
                fragMeans[i] = Double.parseDouble(mu[i]);
            }
        } else {
            fragMeans[0] = 50.0;
            fragMeans[1] = 200.0;
            fragMeans[2] = 400.0;
            fragMeans[3] = 600.0;
        }
        if (stddevs != null) {
            String[] std = stddevs.split(",");
            for (int i = 0; i < std.length; i++) {
                fragStddevs[i] = Double.parseDouble(std[i]);
            }
        } else {
            for (int i = 0; i < fragStddevs.length; i++) {
                fragStddevs[i] = 20.0;
            }
        }

        /*
         * Pull the lengths from the read data. Fragment distribution uses fragments with length > 100 to train the 
	 * nucleosome distributions. the short distribution is set to 50 and remains unchanged
	 * Only occurs if EM training occurs, else use the default settings
         */
        if (fragEM) {
            logger.info("Start to run EM...");
            tmpTime = System.currentTimeMillis();

            pullLargeLengths puller = new pullLargeLengths(bam, index, minMapQ, genomeStats, java.util.Arrays.copyOf(fragMeans, 4));
            double[] lengths = puller.getSampledLengths(10);
            double[] weights = puller.getWeights();
            puller = null;

            //Perform EM training
            HMMR_EM em = new HMMR_EM(weights, java.util.Arrays.copyOfRange(fragMeans, 1, 4),
                    java.util.Arrays.copyOfRange(fragStddevs, 1, 4), lengths);
            em.learn();
            double[] tempMeans = em.getMeans();
            double[] tempLam = em.getLamda();
            em = null;
            for (int i = 0; i < tempMeans.length; i++) {
                //This will update the parameters IFF they were updated. If they become NaN, leave as default
                if (!Double.isNaN(tempMeans[i]) && !Double.isNaN(tempLam[i])) {
                    fragMeans[i + 1] = tempMeans[i];
                    fragStddevs[i + 1] = tempLam[i];
                }
            }
            em = null;
            lengths = null;
            tempMeans = null;
            tempLam = null;
            logger.log(Level.INFO, "EM finished in {0}s", (System.currentTimeMillis() - tmpTime) / 1000);
        }

        logger.info("Fragment length distribution from EM:");
        logger.log(Level.INFO, "NFR: Mean\t{0}\tStdDevs\t{1}", new Object[]{fragMeans[0], fragStddevs[0]}) ;
        logger.log(Level.INFO, "1Ns: Mean\t{0}\tStdDevs\t{1}", new Object[]{fragMeans[1], fragStddevs[1]}) ;
        logger.log(Level.INFO, "NFR: Mean\t{0}\tStdDevs\t{1}", new Object[]{fragMeans[2], fragStddevs[2]}) ;
        logger.log(Level.INFO, "NFR: Mean\t{0}\tStdDevs\t{1}", new Object[]{fragMeans[3], fragStddevs[3]}) ;

        /*
		 * Generate genome wide pileup of read coverage and simultaneously calculate mean and standard deviation
		 * to calculate fold change and zscore. Convert pileup to bedgraph for easier storage. Calculate the pileup
		 * in chunks for memory efficiency
		 * 
		 * NOTE: this method does not currently exist. Instead we use an inputed big wig file to find fold change 
		 * training regions and zscored masked regions
         */
 /*
		 * Below is method to use BigWig input file
		 * Find regions with fold change within determined range to use as training sites.
		 * Find regions with zscore values above certain cutoff to exclude from viterbi.
         */
        Hmm<ObservationVector> hmm = null;

        // Comment out 10_13_18 to eliminate bigwig requirement
        //FoldChange fc = new FoldChange(bigWig,upper,lower,genomeStats);
        //PullWigAboveCutoff z = new PullWigAboveCutoff(bigWig,zscore,genomeStats);
        logger.info("Start to run pileup...");
        tmpTime = System.currentTimeMillis();

        pileup pileupData = new pileup(new SplitBed(genomeStats, vitWindow).getResult(), 0, bam, index, 0);
//		pileup pileupData = new pileup(genomeStats, 0, bam, index, minMapQ);
        bedGraphMath fc = new bedGraphMath(pileupData.getBedGraph());
        logger.info("Building pileups finished in " + (System.currentTimeMillis() - tmpTime) / 1000 + "s");

        pileupData = null;

        //ArrayList<TagNode> bdg2 = fc.getBedGraph();
        //for (int i = 0; i < bdg2.size();i++){
        //	System.out.println(bdg2.get(i).toString3()+"\tBedgraph Signal");
        //}
        double genomeMean = fc.getMean();
        double genomeStd = fc.getSTD();

        ArrayList<TagNode> train = new MergeBed(fc.getBetweenRanges(upper, lower)).getResults();

        ArrayList<TagNode> newTrain = new ArrayList<>();
        int maxTrain;
        if (train.size() > 1000) {
            maxTrain = 1000;
        } else {
            maxTrain = train.size();
        }

        for (int i = 0; i < maxTrain; i++) {
            newTrain.add(train.get(i));
        }

        train = newTrain;
        train = new ExtendBed(train, 5000).getResults();
        //fc = null; 10_13_18
        ArrayList<TagNode> exclude = new MergeBed(fc.getAboveZscore(zscore)).getResults();
        ArrayList<TagNode> addBack = exclude;
        //z=null; 10_13_18

        if (blacklist != null) {
            exclude.addAll(black);
            exclude = new MergeBed(exclude).getResults();//added 8-3-18
        }

        //for (int c = 0; c < exclude.size(); c++) {
        //    logger.info(exclude.get(c).toString() + "\t" + "exclude");
        //}
        newTrain = new ArrayList<TagNode>();
        for (int i = 0; i < train.size(); i++) {
            int counter = 0;
            TagNode node1 = train.get(i);
            for (int a = 0; a < exclude.size(); a++) {
                TagNode node2 = exclude.get(a);
                if (SubtractBed.overlap(node1, node2).hasHit()) {
                    counter++;
                }
            }
            if (counter == 0) {
                newTrain.add(node1);

            }
        }

        train = newTrain;

        //Below added on 12/13/16
        // Allows user to use training set for model generation
        if (trainingRegions != null) {
            bedFileReader trainReader = new bedFileReader(trainingRegions);
            train = trainReader.getData();
        }
        //Above added 12/13/16

        logger.info("Training Regions found and Zscore regions for exclusion found");

        /*
		 * Create the fragment pileup tracks using the training set and the fragment distribution parameters
         */
        if (modelFile == null) {
            logger.info("Start to train model using BW algorithm...");
            tmpTime = System.currentTimeMillis();

            logger.info("Training on " + train.size() + " regions");

            //for (int i = 0; i < train.size(); i++) {
            //    logger.info(train.get(i).toString() + "\t" + "training");
            //}
            FragPileupGen gen = new FragPileupGen(bam, index, train, mode, fragMeans, fragStddevs, minMapQ);
            TrackHolder holder = new TrackHolder((gen.transformTracks(gen.getAverageTracks())), trim);// 8/30/16 transformed tracks 
            //	7/16/18 transformation removed after testing showed it has little effect with new weighted procedure 

            gen = null;

            logger.info("Training Fragment Pileup completed");

            /*
		 * Create the initial model using KMeans and then refine it using Baum-Welch
             */
            KMeansToHMM kmeans = new KMeansToHMM(holder.getDataSet(), k, Integer.MAX_VALUE, true, true, true);
            logger.info("Kmeans Model:\n" + kmeans.getHMM().toString()); // added 7-13-18
            //System.out.println(kmeans.getHMM().toString());

            hmm = new BaumWelch(kmeans.getHMM(), holder.getBWObs(), 150, 0.001).build();

            //System.out.println(hmm.toString());
            kmeans = null;
            holder = null;
            logger.info("Learn HMM with BaumWelch finished in " + (System.currentTimeMillis() - tmpTime) / 1000 + "s");

        } else {
            /*
		 * Use input model if available
             */

            hmm = (Hmm<ObservationVector>) HmmBinaryReader.read(new FileInputStream(modelFile));
        }

        /*
		 * Identify peak state as the state with the highest short read signal.
		 * Identify flanking nucleosome state as the state with the second highest mono read signal.
         */
        logger.info("Start to refine models...");
        tmpTime = System.currentTimeMillis();
        int peak = -1;
        double max = 0.0;
        for (int i = 0; i < hmm.nbStates(); i++) {
            OpdfMultiGaussian pdf = (OpdfMultiGaussian) hmm.getOpdf(i);
            double sh = pdf.mean()[0];
            if (sh > max) {
                peak = i;
                max = sh;
            }
        }
        max = 0.0;
        for (int i = 0; i < hmm.nbStates(); i++) {
            if (i != peak) {
                OpdfMultiGaussian pdf = (OpdfMultiGaussian) hmm.getOpdf(i);
                double mono = pdf.mean()[1];
                if (mono > max) {
                    max = mono;
                }
            }
        }

        /*
		 * Output binary model file
         */
        File outputModel = new File(output + ".model");
        FileOutputStream outModel = new FileOutputStream(outputModel);
        HmmBinaryWriter.write(outModel, hmm);
        outModel.close();
        logger.info("Model created and refined. See " + output + ".model");
        logger.info("Model:\n" + hmm.toString());
        logger.info("Refining model finished in " + (System.currentTimeMillis() - tmpTime) / 1000 + "s");

        /*
		 * Stop program if only model is desired
         */
        if (stopAfterModel) {
            System.exit(0);
        }

        /*
		 * Split the genome file into smaller 25MB chunks 
		 * Can also split into whatever sized chunks the users prefers
		 * May be necessary to split into smaller chunks for machines with less memory
         */
        ArrayList<TagNode> split = new SplitBed(genomeStats, vitWindow).getResult();
        genomeStats = null;

        /*
		 * Subtract excluded regions from the split genome for Viterbi
         */
        ArrayList<TagNode> vitBed = new SubtractBed(split, exclude).getResults();
        split = null;
        exclude = null;

        logger.info("Genome split and subtracted masked regions");

        /*
		 * Run viterbi on the whole genome
         */
        //PrintStream out = new PrintStream(output+".pileup");
        logger.info("Start to decode genome using Viterbi...");
        tmpTime = System.currentTimeMillis();
        ArrayList<TagNode> genomeAnnotation = new ArrayList<TagNode>();
        for (int i = 0; i < vitBed.size(); i++) {
            if (vitBed.get(i).getLength() >= 10) {
                ArrayList<TagNode> tempBed = new ArrayList<TagNode>();
                tempBed.add(vitBed.get(i));
                FragPileupGen vGen = new FragPileupGen(bam, index, tempBed, mode, fragMeans, fragStddevs, minMapQ);
                TrackHolder vHolder = new TrackHolder(vGen.transformTracks(vGen.getAverageTracks()), trim);// 8/30/16 transformed tracks

                //System.out.println(vitBed.get(i).toString()+"\t"+"viterbi bed");
                vGen = null;

                RobustHMM HMM = new RobustHMM(vHolder.getObs(), null, hmm, false, 0, "Vector", 0);
                int[] states = HMM.getStates();
                //Reverse the states for creating peaks
                //ArrayUtils.reverse(states);

                int start = vitBed.get(i).getStart();
                int remainder = vitBed.get(i).getLength() % 10;
                ArrayList<PileupNode2> pile = new ArrayList<PileupNode2>();
                int a;
                for (a = 0; a < states.length - 1; a++) {
                    //PileupNode2 pNode = new PileupNode2(start+(a*10),(double)states[a],vitBed.get(i).getChrom());
                    PileupNode2 pNode = new PileupNode2(start + (a * 10), (double) states[a], vitBed.get(i).getChrom());
                    pile.add(pNode);

                    //out.println(pNode.getChrom()+"\t"+pNode.getBase()+"\t"+(pNode.getBase()+10)+"\t"+"E"+(int)pNode.getScore());
                }
                PileupNode2 pNode = new PileupNode2(start + (((a) * 10) - remainder), (double) states[a], vitBed.get(i).getChrom());
                pile.add(pNode);
                genomeAnnotation.addAll(new PileupToBedGraph(pile, 10).getBedGraph());

                //logger.info(i + " round viterbi done");
            }
        }
        logger.info("Decoding using Viterbi finished in " + (System.currentTimeMillis() - tmpTime) / 1000 + "s");

        //out.close();
        /**
         * Report the final results as peaks, bedgraphs and summits, if desired
         */
        PrintStream bedgraph = null;
        if (bg) {
            bedgraph = new PrintStream(output + ".bedgraph");
        }
        PrintStream pks = null;
        PrintStream summits = null;
        if (peaks) {
            pks = new PrintStream(output + "_peaks.gappedPeak");
            summits = new PrintStream(output + "_summits.bed");
        }
        HashMap<String, ArrayList<TagNode>> bdg = fc.getMappedBedgraph();
//		HashMap<String,ArrayList<TagNode>> bdg = null;

        fc = null;
        HashMap<String, ArrayList<TagNode>> hmmrBdg = bedGraphMath.toMap(genomeAnnotation);
        int counter = 1;
        for (String chr : hmmrBdg.keySet()) {
            ArrayList<TagNode> hmmr = hmmrBdg.get(chr);
            ArrayList<TagNode> signal = bdg.get(chr);
            Collections.sort(hmmr, TagNode.basepairComparator);
            if (signal != null) {
                Collections.sort(signal, TagNode.basepairComparator);
            }
            int index = 0;
            for (int i = 0; i < hmmr.size(); i++) {
                TagNode temp = hmmr.get(i);

                /**
                 * Execute the scoring commands if the state is a peak or if
                 * bedgraph scoring is on
                 */
                if ((int) temp.getScore2() == peak || BGScore) {
                    boolean hasHadOverlap = false;
                    ArrayList<TagNode> overlaps = new ArrayList<TagNode>();
                    for (int a = index; a < signal.size(); a++) {
                        if (SubtractBed.overlap(temp, signal.get(a)).hasHit()) {
                            overlaps.add(signal.get(a));
                            hasHadOverlap = true;
                        } else {
                            if (hasHadOverlap) {
                                index = a;
                                break;
                            }
                        }

                    }
                    ScoreNode scores = bedGraphMath.set(temp, overlaps);
                    if (scoreSys.equals("ave")) {
                        temp.setScore3(scores.getMean());
                    } else if (scoreSys.equals("fc")) {
                        temp.setScore3(scores.getMean() / genomeMean);
                    } else if (scoreSys.equals("zscore")) {
                        temp.setScore3((scores.getMean() - genomeMean) / genomeStd);
                    } else if (scoreSys.equals("med")) {
                        temp.setScore3(scores.getMedian());
                    } else {
                        temp.setScore3(scores.getMax());
                    }
                    if ((int) temp.getScore2() == peak) {
                        temp = bedGraphMath.setSmooth(20, temp, overlaps);
                        temp.setID("Peak_" + counter);
                        if (i > 0) {
                            temp.setUpstream(hmmr.get(i - 1));
                        } else {
                            temp.setUpstream(hmmr.get(i));
                        }
                        if (i < hmmr.size() - 1) {
                            temp.setDownstream(hmmr.get(i + 1));
                        } else {
                            temp.setDownstream(hmmr.get(i));
                        }
                        counter++;
                    }

                }
                /**
                 * report the bedgraph, is desired
                 */
                if (bg) {
                    if (!BGScore) {
                        bedgraph.println(temp.toString2());
                    } else {
                        bedgraph.println(temp.toString_ScoredBdg());
                    }
                }
                /**
                 * report the peaks and summits, if desired
                 */
                if (peaks && (int) temp.getScore2() == peak && temp.getLength() >= minLength) {
                    if (temp.getSummit() != null) {
                        summits.println(temp.getSummit().toString_ScoredSummit());
                    }
                    pks.println(temp.toString_gappedPeak());
                }

            }
        }
        if (bg) {
            bedgraph.close();
        }

        counter = 0;

        for (int i = 0; i < addBack.size(); i++) {
            String chrom = addBack.get(i).getChrom();
            int start = addBack.get(i).getStart();
            int stop = addBack.get(i).getStop();

            pks.println(chrom + "\t" + start + "\t" + stop + "\t" + "HighCoveragePeak_" + counter + "\t.\t.\t0\t0\t255,0,0\t1\t"
                    + addBack.get(i).getLength() + "\t0\t-1\t-1\t-1");
        }
        if (peaks) {
            pks.close();
            summits.close();
        }
        
        Long endTime = System.currentTimeMillis();
        Long total = (endTime - startTime) / 1000;
        logger.info("Total time (seconds)= \t" + total);
    }//main

}
