package projects.pwmbench;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.stream.Stream;

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
		//use AUC-ROC in log'ed intensities
		ROCLOG,
		//use AUC-PR in log'ed intensities
		PRLOG
	}

	private static class Entry{

		private DataSet data;
		private double[] vals;
		private String file;


		public Entry(DataSet data, double[] vals, String file) {
			super();
			this.data = data;
			this.vals = vals;
			this.file = file;
		}



	}


	private static void parseMeme(String file, LinkedList<PFMWrapperTrainSM> models) throws NumberFormatException, IOException, CloneNotSupportedException{
		BufferedReader reader = new BufferedReader( new FileReader(file) );

		String str = null;

		LinkedList<double[]> pwm = new LinkedList<>();
		String name = null;
		boolean inPWM = false;

		while( (str = reader.readLine()) != null ){
			if(str.startsWith("MOTIF ")){
				if(name != null){
					PFMWrapperTrainSM temp = new PFMWrapperTrainSM(DNAAlphabetContainer.SINGLETON, name, pwm.toArray(new double[0][]), 4E-4);
					models.add(temp);
				}
				name = str.replaceAll("^MOTIF ", "");
				pwm.clear();
			}else if(str.startsWith("letter-probability")){
				inPWM = true;
			}else if(str.startsWith("URL") || str.trim().length()==0){
				inPWM = false;
			}else if(inPWM){
				String[] parts = str.trim().split("\\s+");
				double[] line = new double[parts.length];
				if(line.length != 4){
					System.err.println(str);
					throw new RuntimeException("Matrix rows should contain exactly 4 columns");
				}
				for(int i=0;i<line.length;i++){
					line[i] = Double.parseDouble(parts[i]);
				}
				pwm.add(line);
			}
		}

		if(name != null){
			PFMWrapperTrainSM temp = new PFMWrapperTrainSM(DNAAlphabetContainer.SINGLETON, name, pwm.toArray(new double[0][]), 1E-4);
			models.add(temp);

		}

		reader.close();
	}



	private static void parsePhilipp(String file, LinkedList<PFMWrapperTrainSM> models) throws NumberFormatException, IOException, CloneNotSupportedException{
		BufferedReader reader = new BufferedReader( new FileReader(file) );

		String str = null;

		LinkedList<double[]> pwm = new LinkedList<>();
		String name = null;

		while( (str = reader.readLine()) != null ){
			if(str.startsWith(">")){
				if(name != null){
					PFMWrapperTrainSM temp = new PFMWrapperTrainSM(DNAAlphabetContainer.SINGLETON, name, pwm.toArray(new double[0][]), 4E-4);
					models.add(temp);
				}
				name = str.replaceAll(">letter-probability matrix ", "").replaceAll(":.*", "");
				pwm.clear();
			}else{
				String[] parts = str.trim().split("\\s+");
				double[] line = new double[parts.length];
				if(line.length != 4){
					System.err.println(str);
					throw new RuntimeException("Matrix rows should contain exactly 4 columns");
				}
				for(int i=0;i<line.length;i++){
					line[i] = Double.parseDouble(parts[i]);
				}
				pwm.add(line);
			}
		}

		if(name != null){
			PFMWrapperTrainSM temp = new PFMWrapperTrainSM(DNAAlphabetContainer.SINGLETON, name, pwm.toArray(new double[0][]), 4E-4);
			models.add(temp);

		}

		reader.close();
	}



	private static LinkedList<PFMWrapperTrainSM> readModels(String... dFiles ) throws IOException, CloneNotSupportedException{

		LinkedList<PFMWrapperTrainSM> models = new LinkedList<>();

		for(String file:dFiles){
			//parseMeme(file, models);
			parsePhilipp(file, models);

		}



		return models;

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

	private static LinkedList<Entry> readData(String dataDirOrFile) throws IOException, IllegalArgumentException, WrongAlphabetException, EmptyDataSetException{

		Iterator<File> fit = null;
		LinkedList<Entry> all = new LinkedList<>();
		String dataDir = null;

		if( (new File(dataDirOrFile)).isDirectory() ) {
			dataDir = dataDirOrFile;

			Stream<File> stream = Files.walk(Paths.get(dataDirOrFile)).filter(Files::isRegularFile).map(Path::toFile);

			fit = stream.iterator();

		}else {
			File dataFile = new File(dataDirOrFile);
			dataDir = dataFile.getParent();
			if(dataDir == null) {
				dataDir = "";
			}

			LinkedList<File> li = new LinkedList<>();
			li.add(dataFile);

			fit = li.iterator();
		}

		while( fit.hasNext() ){
			File f = fit.next();
			if(f.getAbsolutePath().endsWith(".txt")){

				BufferedReader reader = new BufferedReader(new FileReader(f));

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
				all.add( new Entry( new DataSet("",seqs) , vals.toArray(), f.getAbsolutePath().replaceAll("^"+dataDir+"/?", "")) );

				reader.close();

			}
		}

		return all;

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
				}else {
					neg.add(preds[i]);
				}
			}

			pos.sort();
			neg.sort();

			if(scoring == Scoring.ROC || scoring == Scoring.ROCLOG) {
				AucROC roc = new AucROC();
				return (double)roc.compute(pos.toArray(), neg.toArray()).getResultAt(0).getValue();
			}else {
				AucPR pr = new AucPR();
				return (double)pr.compute(pos.toArray(), neg.toArray()).getResultAt(0).getValue();
			}


		}else {
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
	}


	public static void main(String[] args) throws Exception {

		Scoring scoring = Scoring.valueOf(args[0]);

		LinkedList<Entry> data = readData(args[1]);

		int off = 2;

		String[] sub = new String[args.length-off];
		System.arraycopy(args, off, sub, 0, sub.length);

		LinkedList<PFMWrapperTrainSM> models = readModels(sub);


		Iterator<PFMWrapperTrainSM> modIt = models.iterator();

		Iterator<Entry> datIt = data.iterator();

		while(datIt.hasNext()){
			Entry en = datIt.next();
			System.out.print("\t"+en.file);
		}
		System.out.println();
		while(modIt.hasNext()){
			PFMWrapperTrainSM model = modIt.next();

			System.out.print(model.getName());

			datIt = data.iterator();

			while(datIt.hasNext()){
				Entry en = datIt.next();
				double score = score(model, en.data, en.vals, scoring);
				System.out.print("\t"+score);
			}
			System.out.println();
		}
	}
}
