package buildTreesFromSequenceData;

import java.io.File;

import util.DNAmanipulations;
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
					char twobaseinv=DNAmanipulations.reverse(twobase+"").charAt(0);
					if(onebase!=twobase&&onebase!=twobaseinv){
						System.out.println("Problem with pos: "+i);
						System.out.println(one.getColumn(i)+"\n"+two.getColumn(i));
					}
				}
			}
		}
	}
}
