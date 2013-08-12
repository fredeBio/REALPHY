package buildTreesFromSequenceData;

import java.io.File;

import util.Fasta;
import util.phylogenetics.Alignment;

public class CompareAlignments {
	public static void main(String args[]){
		Alignment one=new Alignment(Fasta.readFasta(new File(args[0])));
		Alignment two=new Alignment(Fasta.readFasta(new File(args[1])));
		for(int i=0;i<one.getLength();i++){
			for(int j=0;j<one.getNumSeqs();j++){
				char onebase=one.getBase(j, i);
				if(onebase!='-'){
					char twobase=two.getBase(j, i);
					if(onebase!=twobase){
						System.out.println("Problem with pos: "+i);
					}
				}
			}
		}
	}
}
