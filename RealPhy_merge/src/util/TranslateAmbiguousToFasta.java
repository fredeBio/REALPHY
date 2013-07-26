package util;

import java.io.*;
import java.util.ArrayList;


public class TranslateAmbiguousToFasta {
	public static void main(String[] args){
		File in=new File(args[0]);
		File out=new File(args[1]);
		ArrayList<Fasta> fasta=Fasta.readFasta(in);
		Fasta.write(translate(fasta), out);
	}
	
	public static ArrayList<Fasta> translate(ArrayList<Fasta> in){
		ArrayList<Fasta> translated=new ArrayList<Fasta>();
		for(int i=0;i<in.size();i++){
			ArrayList<StringBuilder> seqs=translateSeq(in.get(i).getSequence());
			for(int j=0;j<seqs.size();j++){
				translated.add(new Fasta(in.get(i).ident+"_"+j,seqs.get(j).toString()));
			}
		}
		return translated;
	}
	
	public static ArrayList<StringBuilder> translateSeq(String seq){
		ArrayList<StringBuilder> results=new ArrayList<StringBuilder>();
		results.add(new StringBuilder());
		for(int i=0;i<seq.length();i++){
			if(seq.charAt(i)=='['){
				ArrayList<Character> list=new ArrayList<Character>();
				i++;
				while(seq.charAt(i)!=']'){
					list.add(seq.charAt(i));
					i++;
				}
				results=appendChars(results,list);
			}else{
				for(int j=0;j<results.size();j++){
					results.get(j).append(seq.charAt(i));
				}
			}
		}
		return results;
	}
	
	public static ArrayList<StringBuilder> appendChars(ArrayList<StringBuilder> old,ArrayList<Character> list){
		ArrayList<StringBuilder> newList=new ArrayList<StringBuilder>();
		for(int i=0;i<list.size();i++){
			for(int j=0;j<old.size();j++){
				StringBuilder temp=new StringBuilder(old.get(j));
				newList.add(temp.append(list.get(i)));
			}
		}
		return newList;
	}
	
}
