package util;

import java.io.File;
import java.util.ArrayList;
import java.util.BitSet;
public class Find {
	public static void main(String[] args){
		String word=args[0];
		File genome=new File(args[1]);
		StringBuilder[] sb=ReadFasta.readFasta(genome).values().toArray(new StringBuilder[0]);
		int sum=0;
		for(int i=0;i<sb.length;i++){
			String sequence=sb[i].toString();
			sum+=getOccurrences(word,sequence.toUpperCase());
		}
		System.out.println(sum);
		
	}
	public static ArrayList<Integer> getPositions(ArrayList<String> words,ArrayList<String> sequence){
		ArrayList<Integer> pos=new ArrayList<Integer>();
		System.out.println(words.size()+" words.");
		for(int i=0;i<words.size();i++){
			String word=words.get(i);
			for(int j=0;j<sequence.size();j++){
				pos.addAll(getPositions(word, sequence.get(j)));
				pos.addAll(getPositions(DNAmanipulations.reverse(word), sequence.get(j)));

			}
		}
		return pos;
	}
	public static ArrayList<Integer> getPositions(String word,String sequence){
		ArrayList<Integer> pos=new ArrayList<Integer>();
		int i=0;
		while((i=sequence.indexOf(word, i))!=-1){
			pos.add(i);
			i++;
		}

		return pos;
	}

	public static int getOccurrences(ArrayList<BitSet> words,ArrayList<String> sequence){
		int sum=0;
		System.out.println(words.size()+" words.");
		for(int i=0;i<words.size();i++){
			String word=DNAmanipulations.decodeDNA(words.get(i));
			for(int j=0;j<sequence.size();j++){
				sum+=getOccurrences(word, sequence.get(j));
			}
		}
		return sum;
	}
//	public static int getOccurrencesBitSets(ArrayList<BitSet> words,ArrayList<BitSet> sequence){
//		int sum=0;
//		System.out.println(words.size()+"words.");
//		for(int i=0;i<words.size();i++){
//			for(int j=0;j<sequence.size();j++){
//				sum+=getOccurrences(words.get(i), sequence.get(j));
//			}
//		}
//		return sum;
//	}
	public static int getOccurrences(String word,String sequence){
		int i=0;
		int from=0;
		sequence=sequence.toUpperCase();
		word=word.toUpperCase();
		while((from=sequence.indexOf(word,from))!=-1){
			i++;
			from++;
		}
		from=0;
		String rev=DNAmanipulations.reverse(word);
		while((from=sequence.indexOf(rev,from))!=-1){
			i++;
			from++;
		}
		
		return i;
	}
//	public static int getOccurrences(BitSet word,BitSet sequence){
//		int occ=0;
//		if(word.length()<=sequence.length()){
//			for(int i=0;i<sequence.length()-word.length()+1;i++){
//				
//				if(sequence.get(i, i+word.length())==word){
//					occ++;
//				}
//			}
//		}
//		
//		return occ;
//	}
//	public static ArrayList<Integer> getPositions(BitSet word,BitSet sequence){
//		ArrayList<Integer> pos=new ArrayList<Integer>();
//		if(word.size()<=sequence.size()){
//			for(int i=0;i<sequence.size()-word.size();i+=2){
//				BitSet w1=sequence.get(i, i+word.size()-1);
//				w1.set(w1.size());
//				System.out.println(DNAmanipulations.decodeDNA(word));
//				System.out.println(DNAmanipulations.decodeDNA(w1));
//				
//				if(sequence.get(i, i+word.size()-1).equals(word.get(0, word.size()-1))){
//					pos.add(i);
//				}
//			}
//		}
//		
//		return pos;
//	}
//	public static ArrayList<Integer> getPositions(ArrayList<BitSet> words,ArrayList<BitSet> sequence){
//		ArrayList<Integer> pos=new ArrayList<Integer>();
//		System.out.println(words.size()+"words.");
//		for(int i=0;i<words.size();i++){
//			for(int j=0;j<sequence.size();j++){
//				pos.addAll(getPositions(words.get(i), sequence.get(j)));
//				pos.addAll(getPositions(DNAmanipulations.reverse(words.get(i)),sequence.get(j)));
//			}
//		}
//		return pos;
//	}
	
	
	
}
