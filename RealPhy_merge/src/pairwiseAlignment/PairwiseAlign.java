package pairwiseAlignment;

import java.io.*;
import java.util.*;

import statistics.Stats;
import util.*;

public class PairwiseAlign {
	public static void main(String args[]){
		File f1=new File(args[0]);
		File matrix=new File(args[1]);
		double GapOpen=Double.parseDouble(args[2]);
		double GapC=Double.parseDouble(args[3]);
		ArrayList<Fasta> fas1=Fasta.readFasta(f1);
		ArrayList<Double> pw=new ArrayList<Double>();
		HashMap<Character,HashMap<Character,Integer>> subMat=NeedlemanWunsch.readSimilarityMatrix(matrix);
		for(int i=0;i<fas1.size();i++){
			for(int j=i+1;j<fas1.size();j++){
				String seq1=fas1.get(i).getSequence().toUpperCase();
				String seq2=fas1.get(j).getSequence().toUpperCase();
				NeedlemanWunsch nw=new NeedlemanWunsch(seq1,seq2,subMat,GapOpen,GapC);
				double pw1=nw.getPairwiseIdentity();
				nw=new NeedlemanWunsch(DNAmanipulations.reverse(seq1).toUpperCase(),seq2,subMat,GapOpen,GapC); 
				double pw2=nw.getPairwiseIdentity();
				System.out.println(pw1+" "+pw2);
				double identity=Math.max(pw1,pw2);
				pw.add(identity);
				System.out.println(identity);
				System.out.println(nw.getAlignments());
			}
		}
		Stats stats=new Stats(pw);
		
		System.out.println("Average: "+stats.getAverage());
		System.out.println("Standard deviation: "+stats.getStandardDeviation());
		System.out.println("Standard error: "+stats.getStandardError());
		
	}

	
	
}
