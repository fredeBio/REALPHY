package buildTreesFromSequenceData;

import java.io.*;
import java.util.*;

import util.Fasta;

public class MoveSeqToTop {
	public static void main(String args[]){
		File fas=new File(args[0]);
		String id=args[1];
		String name=fas.getName().split("\\.")[0];
		File out=new File(fas.getParent()+"/"+name+"_"+id+".fas");
		move(id,fas,out);
	}
	
	public static int move(String id,File fas,File out){
		ArrayList<Fasta> fasta=Fasta.readFasta(fas);
		if(id.length()==0&&fasta.size()>0){
			Fasta.write(fasta,out );
			return fasta.get(0).getSequence().length();
		}else if(fasta.size()==0){
			return -1;
		}
		for(int i=0;i<fasta.size();i++){
			
			if(fasta.get(i).getIdent().equals(id)){
				Fasta temp=new Fasta(fasta.get(i).getIdent(),fasta.get(i).getSequence());
				fasta.set(i,fasta.get(0));
				fasta.set(0,temp);
			}
			
		}
		Fasta.write(fasta,out );
		return fasta.get(0).getSequence().length();
		
	}
}
