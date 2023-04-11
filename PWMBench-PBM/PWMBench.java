import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

import org.json.simple.JSONObject;

import de.jstacs.classifiers.performanceMeasures.AucPR;
import de.jstacs.classifiers.performanceMeasures.AucROC;
import de.jstacs.data.DataSet;
import de.jstacs.data.EmptyDataSetException;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.alphabets.DNAAlphabetContainer;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.sequenceScores.statisticalModels.trainable.PFMWrapperTrainSM;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.Normalisation;
import de.jstacs.utils.ToolBox;

public class PWMBench {

	public enum Scoring{
		//use intensities and log-sum-occupancy
		ASIS,
		//use intensities and sum-occupancy
		EXP,
		//use log-intensities and log-sum-occupancy
		LOG,
		//use AUC-ROC
		ROC,
		//use AUC-PR
		PR,
		//use AUC-ROC on log'ed intensities
		ROCLOG,
		//use AUC-PR on log'ed intensities
		PRLOG,
		//use Pearson correlation on 8-mers
		MERS,
		//use Pearson correlation on 8-mers for log-intensities
		LOGMERS
	}

	private static int MER = 8;

    private static class Entry{
        private DataSet data;
        private double[] vals;

        public Entry(DataSet data, double[] vals) {
            super();
            this.data = data;
            this.vals = vals;
        }
    }

    private static void pcm_to_pfm(double[][] matrix) {
        for (double[] row : matrix) {
            Normalisation.sumNormalisation(row);
        }
    }

    private static PFMWrapperTrainSM parsePlainMotif(String filename) throws NumberFormatException, IOException, CloneNotSupportedException {
        BufferedReader reader = new BufferedReader( new FileReader(filename) );
        LinkedList<double[]> pcm_list = new LinkedList<>();
        String str = null;
        String name = null;
        while( (str = reader.readLine()) != null ){
            if (str.startsWith(">")) {
                name = str.replaceAll(">", "").trim();
            } else if (str.trim().length()!=0) {
                String[] parts = str.trim().split("\\s+");
                double[] line = new double[parts.length];
                if (line.length != 4) {
                    System.err.println(str);
                    throw new RuntimeException("Matrix rows should contain exactly 4 columns");
                }
                for (int i = 0; i < line.length; i++) {
                    line[i] = Double.parseDouble(parts[i]);
                }
                pcm_list.add(line);
            }
        }
        reader.close();
        double[][] pcm = pcm_list.toArray(new double[0][]);
        pcm_to_pfm(pcm);
        return new PFMWrapperTrainSM(DNAAlphabetContainer.SINGLETON, name, pcm, 4E-4);
    }

    private static int numCol(String str){
        String[] test = str.split("\\s+");
        if(test.length != 2){
            return -1;
        }
        int num = -1;
        try{
            Double.parseDouble(test[1]);
            num = 1;
        }catch(NumberFormatException ex){
            num = -1;
        }
        if(num<0){
            try{
                Double.parseDouble(test[0]);
                num = 0;
            }catch(NumberFormatException ex){
                num = -1;
            }
        }
        return num;
    }

    private static Entry readData(String filename) throws IOException, IllegalArgumentException, WrongAlphabetException, EmptyDataSetException{
        File dataFile = new File(filename);
        if( !dataFile.isFile() ) {
            throw new RuntimeException("Data should be a regular file");
        }
        BufferedReader reader = new BufferedReader(new FileReader(dataFile));

        DoubleList vals = new DoubleList();
        LinkedList<Sequence> seqs = new LinkedList<>();

        String str = reader.readLine();

        int numCol = numCol(str);
        if(numCol<0){
            str = reader.readLine();
            numCol = numCol(str);
        }
        if(numCol<0){
            System.err.println(str);
            throw new RuntimeException("Incorrect data format");
        }

        do{
            if(str.trim().length()>0){
                String[] parts = str.split("\\s+");
                if(parts.length==2){
                    double val = Double.parseDouble(parts[numCol]);
                    Sequence seq = Sequence.create(DNAAlphabetContainer.SINGLETON, parts[1-numCol].trim());
                    vals.add(val);
                    seqs.add(seq);
                }
            }
        }while( (str = reader.readLine()) != null );

        reader.close();

        return new Entry(new DataSet("", seqs), vals.toArray());
    }

    private static double[] getPredictions(PFMWrapperTrainSM model, DataSet data) throws Exception {
        double[] preds = new double[data.getNumberOfElements()];
        for(int i = 0; i< data.getNumberOfElements(); i++){

            Sequence seq = data.getElementAt(i);

            double[] temp = new double[(seq.getLength()- model.getLength()+1)*2];

            for(int j=0,k=0;j<2;j++){
                for(int l = 0; l<seq.getLength()- model.getLength()+1; l++,k++){
                    temp[k] = model.getLogScoreFor(seq, l);
                }
                seq = seq.reverseComplement();
            }

            preds[i] = Normalisation.getLogSum(temp);
        }
        return preds;
    }

    private static double score(PFMWrapperTrainSM model, DataSet data, double[] vals, Scoring scoring) throws Exception {
        double[] preds = getPredictions(model, data);
        return score(model, data, vals, scoring, preds);
    }

    private static double score(PFMWrapperTrainSM model, DataSet data, double[] vals, Scoring scoring, double[] preds) throws Exception {
        if(scoring == Scoring.ROC || scoring == Scoring.PR || scoring == Scoring.ROCLOG || scoring == Scoring.PRLOG) {
            vals = vals.clone();
            if(scoring == Scoring.ROCLOG || scoring == Scoring.PRLOG) {
                double mi = ToolBox.min(vals);
                for(int i=0;i<vals.length;i++) {
                    vals[i] = Math.log(vals[i]-mi+1.0);
                }
            }
            double mean = ToolBox.mean(0, vals.length, vals);
            double sd = ToolBox.sd(0, vals.length, vals);

            double t = mean + 4d * sd;

            double[] vals2 = vals.clone();
            Arrays.sort(vals2);
            double t2 = vals2[vals2.length-50];
            if(t2 < t) {
                t = t2;
            }

            DoubleList pos = new DoubleList();
            DoubleList neg = new DoubleList();

            for(int i=0;i<vals.length;i++) {
                if(vals[i] >= t) {
                    pos.add(preds[i]);
                } else {
                    neg.add(preds[i]);
                }
            }

            pos.sort();
            neg.sort();

            if(scoring == Scoring.ROC || scoring == Scoring.ROCLOG) {
                AucROC roc = new AucROC();
                return (double)roc.compute(pos.toArray(), neg.toArray()).getResultAt(0).getValue();
            } else {
                AucPR pr = new AucPR();
                return (double)pr.compute(pos.toArray(), neg.toArray()).getResultAt(0).getValue();
            }

        } else if(scoring == Scoring.EXP || scoring == Scoring.LOG || scoring == Scoring.ASIS) {
			if(scoring == Scoring.EXP) {
				double mi = ToolBox.min(preds);
				for(int i=0;i<preds.length; i++) {
					preds[i] = Math.exp(preds[i]-mi);
				}
			} else if(scoring == Scoring.LOG) {
                vals = vals.clone();
                double mi = ToolBox.min(vals);
                for (int i = 0; i < vals.length; i++) {
                    vals[i] = Math.log(vals[i] - mi + 1.0);
                }
            }
            return ToolBox.pearsonCorrelation(vals, preds);
		} else {
			HashMap<String, DoubleList[]> merMap = new HashMap<String, DoubleList[]>();
			double mi = ToolBox.min(vals);

			for(int i=0;i<data.getNumberOfElements();i++) {
				Sequence seq = data.getElementAt(i);
				double pred = preds[i];
				double val = (scoring == Scoring.LOGMERS ? Math.log(vals[i]-mi+1.0) : vals[i]);
				for(int j=0;j<seq.getLength()-MER+1;j++) {
					Sequence temp = seq.getSubSequence(j, MER);
					if(temp.toString().compareTo(temp.reverseComplement().toString())>0) {
						temp = temp.reverseComplement();
					}
					String tempStr = temp.toString();
					if(!merMap.containsKey(tempStr)) {
						merMap.put(tempStr, new DoubleList[] {new DoubleList(),new DoubleList()});
					}
					DoubleList[] li = merMap.get(tempStr);
					li[0].add(pred);
					li[1].add(val);
				}
			}

			double[] meanPreds = new double[merMap.size()];
			double[] meanVals = new double[meanPreds.length];

			int i=0;
			for(String key : merMap.keySet()) {
				DoubleList[] li = merMap.get(key);
				meanPreds[i] = li[0].mean(0, li[0].length());
				meanVals[i] = li[1].mean(0, li[1].length());
				i++;
			}

			return ToolBox.pearsonCorrelation(meanVals, meanPreds);
		}
    }

    public static void printScoringMetrics(Entry entry, PFMWrapperTrainSM model, String scoringMode) throws Exception {
        String[] scoringListStrings = scoringMode.split(",");
        if (scoringMode.equals("all") || (scoringListStrings.length > 1)) {
            List<Scoring> scoringList;
            if (scoringMode.equals("all")) {
                scoringList = Arrays.asList(Scoring.values());
            } else {
                scoringList = new ArrayList<>();
                for (String singleScoringMode: scoringListStrings) {
                    scoringList.add(Scoring.valueOf(singleScoringMode));
                }
            }
            Map<String, Double> result = new HashMap<String, Double>();
            double[] preds = getPredictions(model, entry.data);
            for (Scoring scoring : scoringList) {
                double score = score(model, entry.data, entry.vals, scoring, preds.clone());
                result.put(scoring.toString(), score);
            }
            System.out.println(new JSONObject(result));
        } else {
            Scoring scoring = Scoring.valueOf(scoringMode);
            double score = score(model, entry.data, entry.vals, scoring);
            System.out.println(score);
        }
    }

    public static void main(String[] args) throws Exception {
        String datasetFilename = args[1];
        Entry entry = readData(datasetFilename);
        if (args[2].equals("-")) {
            Scanner scanner = new Scanner(System.in);
            while (scanner.hasNextLine()) {
                String motifFilename = scanner.nextLine();
                PFMWrapperTrainSM model = parsePlainMotif(motifFilename);
                System.out.print(datasetFilename + "\t" + motifFilename + "\t");
                printScoringMetrics(entry, model, args[0]);
            }
            scanner.close();
        } else {
            PFMWrapperTrainSM model = parsePlainMotif(args[2]);
            printScoringMetrics(entry, model, args[0]);
        }
    }
}
