package WigMath;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math3.stat.descriptive.rank.Max;

//import HMMR_ATAC.SplitBed;
//import HMMR_ATAC.SubtractBed;
import Node.ScoreNode;
import Node.TagNode;

public class bedGraphMath {

    private HashMap<String, ArrayList<TagNode>> bedgraph;
//	private HashMap<Integer,ArrayList<TagNode>> bedgraph2;
//	private ArrayList<TagNode> genomeStats;
//	private ArrayList<TagNode> index;

    private double mean;
    private double std;
    private double score;
    private double median;
    private double max;
    private ArrayList<Double> val;
    private TagNode node;

    public bedGraphMath(ArrayList<TagNode> bdg) {
        bedgraph = toMap(bdg);

        setMeanAndStd();
    }

//	public bedGraphMath(ArrayList<TagNode> bdg,ArrayList<TagNode> genome){
//		genomeStats = genome;
//		index = makeIndex(25000000);
//		bedgraph2 = toMap2(setIndexArray(index,bdg));
//		bedgraph = toMap(bdg);
//		setMeanAndStd();
//	}
//	
//	private HashMap<Integer,ArrayList<TagNode>> toMap2(ArrayList<TagNode> i){
//		HashMap<Integer,ArrayList<TagNode>> map = new HashMap<Integer,ArrayList<TagNode>>();
//		for (int x = 0;x < i.size();x++){
//			int index = i.get(x).getScore();
//			if (map.containsKey(index)){
//				ArrayList<TagNode> temp = map.get(index);
//				temp.add(i.get(x));
//				map.put(index, temp);
//			}
//			else{
//				ArrayList<TagNode> temp = new ArrayList<TagNode>();
//				temp.add(i.get(x));
//				map.put(index, temp);
//			}
//		}
//		
//		return map;
//		
//	}
//	
//	
//	private ArrayList<TagNode> setIndexArray(ArrayList<TagNode> index,ArrayList<TagNode> input){
//		ArrayList<TagNode> temp = new ArrayList<TagNode>();
//		for (int i = 0;i < input.size();i++){
//			temp.addAll(setIndex(input.get(i),index));
//		}
//		return temp;
//	}
//	private ArrayList<TagNode> setIndex(TagNode t,ArrayList<TagNode> index){
//		ArrayList<TagNode> temp = new ArrayList<TagNode>();
//		for (int a = 0;a < index.size();a++){
//			TagNode node2 = index.get(a);
//			if (SubtractBed.overlap(t, node2).hasHit()){
//				t.setScore(node2.getScore());
//				temp.add(t);
//			}
//		}
//		return temp;
//	}
//	private ArrayList<TagNode> makeIndex(int window){
//		ArrayList<TagNode> index = new ArrayList<TagNode>();
//		index = new SplitBed(genomeStats,window).getResult();
//		Collections.sort(index,TagNode.basepairComparator);
//		for (int i = 0;i < index.size();i++){
//			index.get(i).setScore(i);
//		}
//		
//		return index;
//	}
//	
    public ArrayList<TagNode> getBedGraph() {
        ArrayList<TagNode> results = new ArrayList<TagNode>();
        for (String chr : bedgraph.keySet()) {

            ArrayList<TagNode> inTemp = bedgraph.get(chr);
            Collections.sort(inTemp, TagNode.basepairComparator);
            for (int i = 0; i < inTemp.size(); i++) {

                results.add(inTemp.get(i));

            }
        }
        return results;
    }

    public double getMean() {
        return mean;
    }

    public double getSTD() {
        return std;
    }

    public HashMap<String, ArrayList<TagNode>> getMappedBedgraph() {
        return bedgraph;
    }

    /**
     * Access the average score
     *
     * @return a double representing the average score
     */
    public double getScore() {
        return score;
    }

    /**
     * Access the median score
     *
     * @return a double representing the median score
     */
    public double getMedian() {
        return median;
    }

    /**
     * Access the max score
     *
     * @return a double representing the maximum score
     */
    public double getMax() {
        return max;
    }

    /**
     * Access the bigwig scores
     *
     * @return an ArrayList of doubles representing all the bigwig scores
     */
    public ArrayList<Double> getValues() {
        return val;
    }

    /**
     * Access the position with the maximum score
     *
     * @return a TagNode representing the position having the maximum score
     */
    public TagNode getMaxPos() {
        return node;
    }

    public static TagNode setSmooth(int stdev, TagNode tag, ArrayList<TagNode> overlaps) throws IOException {

        // Use a window size equal to +/- 3 SD's
        double[] filter = new double[6 * stdev + 1];
        double sum = 0;
        for (int i = 0; i < filter.length; i++) {
            double x = i - 3 * stdev;
            double value = (double) Math
                    .exp(-(x * x) / (2 * stdev * stdev));
            filter[i] = value;
            sum += value;
        }
        // Normalize so that the filter is area-preserving (has total area = 1)
        for (int i = 0; i < filter.length; i++) {
            filter[i] /= sum;
        }
        int start = tag.getStart();
        int stop = tag.getStop();
        int paddedStart = start - (3 * stdev);
        int paddedStop = stop + (3 * stdev);
        double[] data = new double[paddedStop - paddedStart];
        Collections.sort(overlaps, TagNode.basepairComparator);
        for (int a = 0; a < overlaps.size(); a++) {
            TagNode node2 = overlaps.get(a);

            double value = node2.getScore2();

            for (int i = node2.getStart(); i < node2.getStop(); i++) {
                if (i >= paddedStart && i < paddedStop) {
                    data[i - paddedStart] = value;
                }
            }

        }
        // Convolve the data with the filter
        double[] smoothed = new double[stop - start];
        for (int i = 0; i < smoothed.length; i++) {
            for (int j = 0; j < filter.length; j++) {
                smoothed[i] += data[i + j] * filter[j];
            }
        }
        double tempMax = 0.0;
        for (int i = 0; i < smoothed.length; i++) {
            if (smoothed[i] > tempMax) {
                tempMax = smoothed[i];
                tag.setSummit(new TagNode(tag.getChrom(), start + i, start
                        + i + 1));
            }
        }

        return tag;

    }

    public static ScoreNode set(TagNode tag, ArrayList<TagNode> overlaps) {

        Max m = new Max();
        Mean mu = new Mean();
        ArrayList<Double> values = new ArrayList<Double>();
        Collections.sort(overlaps, TagNode.basepairComparator);
        for (int a = 0; a < overlaps.size(); a++) {
            TagNode node2 = overlaps.get(a);

            double value = node2.getScore2();
            m.increment(value);
            for (int i = node2.getStart(); i < node2.getStop(); i++) {
                if (i >= tag.getStart() && i < tag.getStop()) {
                    values.add(value);
                    mu.increment(value);

                }
            }

        }
        double ave;
        double median;
        double max;
        if (values.size() > 0) {
            ave = mu.getResult();
            max = m.getResult();
            median = values.get(values.size() / 2);
        } else {
            ave = 0;
            max = 0;
            median = 0;
        }
        ScoreNode output = new ScoreNode(ave, median, max, values);
        return output;

    }

    public ArrayList<TagNode> getAboveZscore(double z) {
        ArrayList<TagNode> results = new ArrayList<TagNode>();
        for (String chr : bedgraph.keySet()) {

            ArrayList<TagNode> inTemp = bedgraph.get(chr);
            Collections.sort(inTemp, TagNode.basepairComparator);
            for (int i = 0; i < inTemp.size(); i++) {
                double value = inTemp.get(i).getScore2();
                if (((value - mean) / std) >= z) {
                    results.add(inTemp.get(i));
                }
            }
        }
        return results;
    }

    public ArrayList<TagNode> getBetweenZRanges(double upperZ, double lowerZ) {
        ArrayList<TagNode> results = new ArrayList<TagNode>();
        for (String chr : bedgraph.keySet()) {
            ArrayList<TagNode> tmpresults = new ArrayList<TagNode>();
            ArrayList<TagNode> inTemp = bedgraph.get(chr);
            Collections.sort(inTemp, TagNode.basepairComparator);
            System.out.println(chr + ": Total number of tags:" + inTemp.size());
            for (int i = 0; i < inTemp.size(); i++) {
                double value = inTemp.get(i).getScore2();
                if (((value - mean) / std) >= lowerZ && ((value - mean) / std) <= upperZ) {
                    tmpresults.add(inTemp.get(i));
                }
            }
            //System.out.println(chr + ": Total number of tags after filtering:" + tmpresults.size());
            if ( tmpresults.size() == 0 ) { continue; }
            // merge overlap
            TagNode first = tmpresults.get(0);
            for (int i = 1; i < tmpresults.size(); i++) {
                TagNode next = tmpresults.get(i);
                if (first.getStop() + 500 >= next.getStart()) {
                    first.setStop(next.getStop());
                } else {
                    results.add(first);
                    first = next;
                }
            }
            results.add(first);
        }
        return results;
    }

    public ArrayList<TagNode> getBetweenRanges(double upper, double lower) {
        ArrayList<TagNode> results = new ArrayList<TagNode>();
        for (String chr : bedgraph.keySet()) {

            ArrayList<TagNode> inTemp = bedgraph.get(chr);
            Collections.sort(inTemp, TagNode.basepairComparator);
            for (int i = 0; i < inTemp.size(); i++) {
                double value = inTemp.get(i).getScore2();
                if ((value / mean) >= lower && (value / mean) <= upper) {
                    results.add(inTemp.get(i));
                }
            }
        }
        return results;
    }

    private void setMeanAndStd() {
        Mean mu = new Mean();
        StandardDeviation dev = new StandardDeviation();
        for (String chr : bedgraph.keySet()) {

            ArrayList<TagNode> inTemp = bedgraph.get(chr);
            Collections.sort(inTemp, TagNode.basepairComparator);
            for (int i = 0; i < inTemp.size(); i++) {
                int length = inTemp.get(i).getLength();
                double value = inTemp.get(i).getScore2();
                for (int a = 0; a < length; a++) {
                    mu.increment(value);
                    dev.increment(value);
                }
            }
        }
        mean = mu.getResult();
        std = dev.getResult();
    }

    public static HashMap<String, ArrayList<TagNode>> toMap(ArrayList<TagNode> i) {
        HashMap<String, ArrayList<TagNode>> map = new HashMap<String, ArrayList<TagNode>>();
        for (int x = 0; x < i.size(); x++) {
            String chr = i.get(x).getChrom();
            if (map.containsKey(chr)) {
                ArrayList<TagNode> temp = map.get(chr);
                temp.add(i.get(x));
                map.put(chr, temp);
            } else {
                ArrayList<TagNode> temp = new ArrayList<TagNode>();
                temp.add(i.get(x));
                map.put(chr, temp);
            }
        }

        return map;

    }
}
