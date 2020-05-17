import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.stream.Stream;

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
		LOG
	}
	
	private static class Entry{

		private DataSet data;
		private double[] vals;
		
		public Entry(DataSet data, double[] vals) {
			super();
			this.data = data;
			this.vals = vals;
		}
	}

	private static PFMWrapperTrainSM parsePlainMotif(String filename) throws NumberFormatException, IOException, CloneNotSupportedException {
		BufferedReader reader = new BufferedReader( new FileReader(filename) );
		LinkedList<double[]> pwm = new LinkedList<>();
		String str = null;
		String name = null;
		while( (str = reader.readLine()) != null ){
			if (str.startsWith(">")) {
				name = str.replaceAll(">", "").trim();
			} else {
				String[] parts = str.trim().split("\\s+");
				double[] line = new double[parts.length];
				if (line.length != 4) {
					System.err.println(str);
					throw new RuntimeException("Matrix rows should contain exactly 4 columns");
				}
				for (int i = 0; i < line.length; i++) {
					line[i] = Double.parseDouble(parts[i]);
				}
				pwm.add(line);
			}
		}
		reader.close();
		return new PFMWrapperTrainSM(DNAAlphabetContainer.SINGLETON, name, pwm.toArray(new double[0][]), 4E-4);
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
	
	private static Entry readData(String dataDirOrFile) throws IOException, IllegalArgumentException, WrongAlphabetException, EmptyDataSetException {
		File dataFile = new File(dataDirOrFile);
		if( !dataFile.isFile() ) {
			throw new RuntimeException("Data should be a regular file");
		}
		BufferedReader reader = new BufferedReader(new FileReader(dataFile));
		
		DoubleList vals = new DoubleList();
		LinkedList<Sequence> seqs = new LinkedList<>();
		
		String str = reader.readLine();
		
		int numCol = numCol(str);
		if (numCol < 0) {
			str = reader.readLine();
			numCol = numCol(str);
		}
		if (numCol < 0) {
			System.err.println(str);
			throw new RuntimeException("Incorrect data format");
		}
		
		do {
			if (str.trim().length() > 0) {
				String[] parts = str.split("\\s+");
				if (parts.length == 2) {
					double val = Double.parseDouble( parts[numCol] );
					Sequence seq = Sequence.create(DNAAlphabetContainer.SINGLETON, parts[1 - numCol].trim());
					vals.add(val);
					seqs.add(seq);
				}
			}
		} while( (str = reader.readLine()) != null );
		
		reader.close();

		return new Entry(new DataSet("", seqs), vals.toArray());
	}
	
	
	private static double score(PFMWrapperTrainSM model, DataSet data, double[] vals, Scoring scoring) throws Exception{
		
		double[] preds = new double[data.getNumberOfElements()];
		for(int i=0;i<data.getNumberOfElements();i++){
			
			Sequence seq = data.getElementAt(i);
			
			double[] temp = new double[(seq.getLength()-model.getLength()+1)*2];
			
			for(int j=0,k=0;j<2;j++){
				for(int l=0;l<seq.getLength()-model.getLength()+1;l++,k++){
					temp[k] = model.getLogScoreFor(seq, l);
				}
				seq = seq.reverseComplement();
			}
			
			preds[i] = Normalisation.getLogSum(temp);
		}
		
		if(scoring == Scoring.EXP) {
			double mi = ToolBox.min(preds);
			for(int i=0;i<preds.length; i++) {
				preds[i] = Math.exp(preds[i]-mi);
			}
		}else if(scoring == Scoring.LOG) {
			vals = vals.clone();
			double mi = ToolBox.min(vals);
			for(int i=0;i<vals.length;i++) {
				vals[i] = Math.log(vals[i]-mi+1.0);
			}
		}
			
		return ToolBox.pearsonCorrelation(vals, preds);
		
	}
	
	
	public static void main(String[] args) throws Exception {
		Scoring scoring = Scoring.valueOf(args[0]);
		Entry entry = readData(args[1]);
		PFMWrapperTrainSM model = parsePlainMotif(args[2]);
		double score = score(model, entry.data, entry.vals, scoring);
		System.out.println(score);
	}

}
