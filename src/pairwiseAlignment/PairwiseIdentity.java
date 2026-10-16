package pairwiseAlignment;

import java.io.*;
import java.util.*;

import statistics.Stats;
import util.Fasta;

public class PairwiseIdentity {
	public static void main(String args[]){
		File folder=new File(args[0]);
		File[] files=folder.listFiles();
		for(int k=0;k<files.length;k++){
			if(!files[k].getAbsolutePath().endsWith("fas")&&!files[k].getAbsolutePath().endsWith("fasta"))continue;
			System.out.println(files[k]);
			ArrayList<Fasta> fas1=Fasta.readFasta(files[k]);
			ArrayList<Double> pw=new ArrayList<Double>();
			try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(files[k]+".out"));
			for(int i=0;i<fas1.size();i++){
				for(int j=i+1;j<fas1.size();j++){
					bw.write(NeedlemanWunsch.getPairwiseIdentity(fas1.get(i).getSequence(), fas1.get(j).getSequence())+"\r\n");
					pw.add(NeedlemanWunsch.getPairwiseIdentity(fas1.get(i).getSequence(), fas1.get(j).getSequence()));
				}
			}
			bw.close();
			}catch(IOException e){
				e.printStackTrace();
				System.exit(-1);
			}
			Stats stats=new Stats(pw);

			System.out.println("Average: "+stats.getAverage());
			System.out.println("Standard deviation: "+stats.getStandardDeviation());
			System.out.println("Standard error: "+stats.getStandardError());
		}
	}
}
