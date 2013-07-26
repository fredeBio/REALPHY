package util;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.Iterator;

public class GenerateMutatedSequences {
	ArrayList<ArrayList<BitSet>> words=new ArrayList<ArrayList<BitSet>>();

	
	public static void main(String args[]){
		File out=new File(args[0]);
		
		int mutations=Integer.parseInt(args[1]);
		ArrayList<BitSet> words=new ArrayList<BitSet>();
		for(int i=2;i<args.length;i++){
			words.add(DNAmanipulations.codeDNA(args[i]));
		}
		GenerateMutatedSequences gm=new GenerateMutatedSequences(words);
		for(int i=0;i<mutations;i++){
			gm.mutateSequences();
		}
		gm.write( out,mutations);
	}
	public GenerateMutatedSequences(ArrayList<BitSet> Words){
		words.add(Words);

	}
	public GenerateMutatedSequences(ArrayList<BitSet> Words,int mut){
		words.add(Words);
		for(int i=1;i<=mut;i++){
			mutateSequences();
		}
		
	}
	public GenerateMutatedSequences(ArrayList<BitSet> Words,int mut,int notMut1,int notMut2){
		words.add(Words);
		for(int i=1;i<=mut;i++){
			mutateSequences(notMut1,notMut2);
		}
		
	}
	public void mutateSequences(){
		ArrayList<BitSet> newList=new ArrayList<BitSet>();
		for(int i=0;i < words.get(words.size()-1).size();i++){
			newList.addAll(mutateSequence(words.get(words.size()-1).get(i)));
		}
		newList.addAll(words.get(words.size()-1));
		words.add(eliminateDoubles(newList));
	}
	public void mutateSequences(int notMut1,int notMut2){
		ArrayList<BitSet> newList=new ArrayList<BitSet>();
		for(int i=0;i < words.get(words.size()-1).size();i++){
			newList.addAll(mutateSequence(words.get(words.size()-1).get(i),notMut1,notMut2));
		}
		newList.addAll(words.get(words.size()-1));
		words.add(eliminateDoubles(newList));
	}
//	private static ArrayList<String> mutateSequence(String seq){
//		char[] word=seq.toCharArray();
//		ArrayList<String> list=new ArrayList<String>();
//		for(int i=0;i<word.length;i++){
//			for(int k=0;k<4;k++){
//				String DNA=DNAmanipulations.decodeDNA(k,1);
//				if(word[i]!=DNA.charAt(0)){
//					word[i]=DNA.charAt(0);
//					list.add(new String(word));
//					word=seq.toCharArray();
//				}
//			}
//		}
//		return list;
//	}
	private static ArrayList<BitSet> mutateSequence(BitSet seq,int notMut1,int notMut2){
		ArrayList<BitSet> list=new ArrayList<BitSet>();
		for(int i=0;i<seq.length()-1;i+=2){
			if(i==notMut1*2||i==notMut2*2)continue;
			for(int k=0;k<4;k++){
				BitSet mut;
				if((mut =DNAmanipulations.mutate(seq, i, k)).length()>0){
					list.add(mut);
				}
				
			}
		}
		return list;
	}
	private static ArrayList<BitSet> mutateSequence(BitSet seq){
		ArrayList<BitSet> list=new ArrayList<BitSet>();
		for(int i=0;i<seq.length()-1;i+=2){
			for(int k=0;k<4;k++){
				BitSet mut;
				if((mut =DNAmanipulations.mutate(seq, i, k)).length()>0){
					list.add(mut);
				}
				
			}
		}
		return list;
	}
	
	public  void write(File out,int mut){
		try{
			BufferedWriter bw = new BufferedWriter(new FileWriter(out));
			for(int i=0;i<words.get(mut).size();i++){
				bw.write(">"+i+"\n"+DNAmanipulations.decodeDNA(words.get(mut).get(i))+"\n");
			}
			bw.close();
		}catch(IOException e){
			System.err.println(e.toString());
		}
	}
	

	
	public ArrayList<BitSet> getList(int mut){
		return words.get(mut);
	}
	
	private  static ArrayList<BitSet> eliminateDoubles(ArrayList<BitSet> words){
		ArrayList<BitSet> newList=new ArrayList<BitSet>();
		HashMap<BitSet,Boolean> temp =new HashMap<BitSet, Boolean>();
		for(int i=0;i<words.size();i++){
			if(!temp.containsKey(DNAmanipulations.reverse(words.get(i))))temp.put(words.get(i), true);
		}
		Iterator<BitSet> it=temp.keySet().iterator();
		while(it.hasNext()){
			newList.add(it.next());
		}
		return newList;
		
	}
	
}
