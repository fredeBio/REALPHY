package buildTreesFromSequenceData;

import java.util.*;

public class Genes {
	String geneID;
	TreeMap<Integer,Character> polymorphisms;
	TreeMap<Integer,Integer> gaps;
	int length;
	int pos;
	char orient;
	String contig;
	public Genes(String GeneID,TreeMap<Integer,Character> Polymorphisms,int Length,int Position,char Orientation,String Contig){
		geneID=GeneID;
		polymorphisms=Polymorphisms;
		length=Length;
		pos=Position;
		orient=Orientation;
		contig=Contig;
		 gaps=new TreeMap<Integer, Integer>();
	}


}
