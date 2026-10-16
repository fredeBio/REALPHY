package pairwiseAlignment;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;


import util.Fasta;
import util.Shuffle;

public class FolderPairwiseIdentity implements FilenameFilter{
		public static void main(String args[]){
			File folder=new File(args[0]);
			File matrix=new File(args[1]);
			int repetitions=Integer.parseInt(args[2]);
			double gapOpen=10;
			double gapCont=0.1;
			File[] files=folder.listFiles(new FolderPairwiseIdentity());
			HashMap<Character,HashMap<Character,Integer>> subMat=NeedlemanWunsch.readSimilarityMatrix(matrix);
			try{
				
			BufferedWriter bw=new BufferedWriter(new FileWriter(new File(folder+"/results.out")));
			BufferedWriter bw2=new BufferedWriter(new FileWriter(new File(folder+"/resultsShuffle.out")));
			bw.write("\t");
			for(int i=0;i<files.length;i++){
				bw.write(files[i].getName()+"\t");
			}
			bw.write("\r\n");
			for(int i=0;i<files.length;i++){
				bw.write(files[i].getName()+"\t");
				for (int m=0;m<i;m++){
					bw.write("\t");
					bw2.write("\t");
				}
				for(int j=i+1;j<files.length;j++){
					if(i==j)continue;
					System.out.print(files[i]+" "+files[j]);
					ArrayList<Fasta> fas1=Fasta.readFasta(files[i]);
					ArrayList<Fasta> fas2=Fasta.readFasta(files[j]);
					if (fas1.size()!=1){
						System.err.println("Only one fasta sequence per file permitted! "+files[i]);
						System.exit(-1);
					}
					if (fas2.size()!=1){
						System.err.println("Only one fasta sequence per file permitted! "+files[j]);
						System.exit(-1);
					}
					String seq1=fas1.get(0).getSequence();
					String seq2=fas2.get(0).getSequence();
					NeedlemanWunsch nw=new NeedlemanWunsch(seq1, seq2,subMat , gapOpen,gapCont);
					double pwi=nw.getPairwiseIdentity();					
					System.out.println(pwi);
					bw.write(pwi+"\t");
					double pValue=getPValue(seq1,seq2,repetitions,subMat,gapOpen,gapCont,pwi);
					bw2.write(pValue+"\t");
					System.out.println(pValue);
				}
				bw.write("\r\n");
				bw2.write("\r\n");

			}
			bw.close();
			bw2.close();
			}catch(IOException e){
				e.printStackTrace();
				System.exit(-1);
			}
			
		}
		
		public static double getPValue(String seq1,String seq2,int reps,HashMap<Character,HashMap<Character,Integer>> subMat,double gapOpen,double gapCont,double threshold){
			int occ=0;
			for(int i=0;i<reps;i++){
				String seq1shuff=Shuffle.shuffle(seq1);
				String seq2shuff=Shuffle.shuffle(seq2);

				NeedlemanWunsch nw=new NeedlemanWunsch(seq1shuff, seq2shuff,subMat , gapOpen,gapCont);
				double pwi=nw.getPairwiseIdentity();
				if(pwi>=threshold){
					occ++;
				}
			}
			return (1.0*occ)/(reps*1.0);
		}

		@Override
		public boolean accept(File dir, String name) {
			// TODO Auto-generated method stub
			if(name.endsWith(".fasta")||name.endsWith(".fas")||name.endsWith(".faa")||name.endsWith(".fna"))return true;
			return false;
		}

		
	
}
