package util;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;



public class DNAmanipulations {
	public static String reverse(String seq){
		StringBuilder sb=new StringBuilder();
		seq=seq.toUpperCase();
		for(int i=seq.length()-1;i>-1;i--){
			sb.append(seq.charAt(i)=='A'?'T':seq.charAt(i)=='T'?'A':seq.charAt(i)=='C'?'G':seq.charAt(i)=='G'?'C':'N');
		}
		
		return sb.toString();
	}
	
	/**
	 * Generates the complement of the sequence. NOT the reverse complement. I.e. it just translates every nucleotide to the complementary nucleotide.
	 * @param seq
	 * @return
	 */
	public static String complement(String seq){
		StringBuilder sb=new StringBuilder();
		seq=seq.toUpperCase();
		for(int i=0;i<seq.length();i++){
			sb.append(seq.charAt(i)=='A'?'T':seq.charAt(i)=='T'?'A':seq.charAt(i)=='C'?'G':seq.charAt(i)=='G'?'C':'N');
		}
		
		return sb.toString();
	}
	
	public static String generateRandomSequence(int length,double GC){
		StringBuffer sb=new StringBuffer();
		for(int i=0;i<length;i++){
			sb.append(getRandomBase(GC));
		}
		return sb.toString();
	}


	
	public static char getRandomBase(double GC){
		double rand=Math.random();
		double r=Math.random();
		if(rand>GC){
			if(r>0.5){
				return 'A';
			}else{
				return 'T';
			}
		}else {
			if(r>0.5){
				return 'C';

			}else{
				return 'G';
			}
		}
	}
	
	public static String randomizeSequenceMutateToSelf(String seq,double mutationProb,double GC){
		StringBuffer sb=new StringBuffer(seq);
		int seqlength=seq.length();
		
		int numberMutations=(int)(seq.length()*mutationProb);
		for(int i=0;i<numberMutations;i++){
			int randPosition=(int)(Math.random()*seqlength);
			sb.setCharAt(randPosition, getMutationSelfMutate(seq.charAt(randPosition),GC));
			
		}
		return sb.toString();
	}
	
	public static String randomizeSequence(String seq,double mutationProb,double GC){
		StringBuffer sb=new StringBuffer(seq);
		int seqlength=seq.length();
		
		int numberMutations=(int)(seq.length()*mutationProb);
		for(int i=0;i<numberMutations;i++){
			int randPosition=(int)(Math.random()*seqlength);
			sb.setCharAt(randPosition, getMutation(seq.charAt(randPosition),GC));
			
		}
		return sb.toString();
	}
	public static String randomizeSequenceSelfMutate(String seq,double mutationProb,double GC){
		StringBuffer sb=new StringBuffer(seq);
		int seqlength=seq.length();
		
		int numberMutations=(int)(seq.length()*mutationProb);
		for(int i=0;i<numberMutations;i++){
			int randPosition=(int)(Math.random()*seqlength);
			sb.setCharAt(randPosition, getMutationSelfMutate(seq.charAt(randPosition),GC));
			
		}
		return sb.toString();
	}
	public static char getMutationSelfMutate(char in,double GC){

			char base=getRandomBase(GC);

				return base;

	}
	
	public static char getMutation(char in,double GC){
		boolean mut=false;
		while(mut!=true){
			char base=getRandomBase(GC);
			if(in!=base){
				return base;
			}
		}
		return 0;
	}
	
	public static BitSet mutate(BitSet seq,int pos,int mutation){
		BitSet mut=new BitSet();
		mut.or(seq);
		if(pos%2==1){
			System.err.println("Position is wrong, has to be dividable by two!");
		}else{
			if(mut.get(pos)!= ((mutation)%2==0) ||mut.get(pos+1)!= ((mutation/2)%2==0)  ){
				mut.set(pos,(mutation)%2==0);
				mut.set(pos+1,(mutation/2)%2==0);
				return mut;
			}
		}
		return new BitSet();
	}
	
	public static BitSet reverse(BitSet code){
		BitSet rev=new BitSet();
		int j=0;
		for(int i=code.length()-3;i>-1;i=i-2){
			if(code.get(i)==false && code.get(i+1)==false){
				//T
				rev.set(j*2,false);
				rev.set(j*2+1,true);
			}else if(code.get(i)==false && code.get(i+1)==true){
				//A
				rev.set(j*2,false);
				rev.set(j*2+1,false);
			}else if(code.get(i)==true && code.get(i+1)==false){
				//G
				rev.set(j*2,true);
				rev.set(j*2+1,true);
			}else if(code.get(i)==true && code.get(i+1)==true){
				//C
				rev.set(j*2,true);
				rev.set(j*2+1,false);
			}
			j++;
			
		}
		rev.set(j*2,true);
		return rev;
	}
	public static BitSet flip(BitSet code){
		BitSet rev=new BitSet();
		int j=0;
		for(int i=code.length()-3;i>-1;i=i-2){
				rev.set(j*2,code.get(i));
				rev.set(j*2+1,code.get(i+1));

			j++;
			
		}
		rev.set(j*2,true);
		return rev;
	}
	public static BitSet append(BitSet b1,BitSet b2){
		int b1size=b1.length();
		int b2size=b2.length();
		if(b1size==0)return b2;

		for(int i=b1size-1;i<b1size+b2size;i++){
			
			b1.set(i,b2.get(i-b1size+1));
		}
		return b1;
	}
	
	public static int next(BitSet reference,BitSet search,int from){
		int pos=-1;
		if(from%2!=0){
			System.err.println("\"From\" has to be dividable by 2!");
			return -1;
		}
		for(int i=from;i<reference.length()-search.length()+1;i+=2){
			if(DNAmanipulations.get(reference,i,i+search.length()-1).equals(search)){
				return i;
			}
		}
		return pos;
	}
	
	public static BitSet get(BitSet b,int start,int end){
		if(end>=b.length())end=b.length()-1;
		if(start>=b.length())return new BitSet();
		BitSet nB=b.get(start,end);
		nB.set(end-start);
		return nB;
	}
/*	public static BitSet complement(BitSet code){
		BitSet rev=new BitSet();
		int j=0;
		for(int i=0;i<code.length()-2;i=i+2){
			if(code.get(i)==false && code.get(i+1)==false){
				//T
				rev.set(j*2,false);
				rev.set(j*2+1,true);
			}else if(code.get(i)==false && code.get(i+1)==true){
				//A
				rev.set(j*2,false);
				rev.set(j*2+1,false);
			}else if(code.get(i)==true && code.get(i+1)==false){
				//G
				rev.set(j*2,true);
				rev.set(j*2+1,true);
			}else if(code.get(i)==true && code.get(i+1)==true){
				//C
				rev.set(j*2,true);
				rev.set(j*2+1,false);
			}
			j++;
			
		}
		rev.set(j*2,true);
		return rev;
	}
	*/
	public  static String translate(String DNA,HashMap<String,String> code){
		StringBuilder AA=new StringBuilder("");
		for(int i=0;i<DNA.length()-2;i+=3){
			AA.append(code.get(DNA.substring(i,i+3).toUpperCase()));
		}
		return AA.toString();
	}
	
	public static HashMap<String,String> code(){
		HashMap<String,String> code=new HashMap<String, String>();
		code.put("TTT", "F");
		code.put("TTC", "F");
		code.put("TTG", "L");
		code.put("TTA", "L");
		code.put("CTT", "L");
		code.put("CTC", "L");
		code.put("CTA", "L");
		code.put("CTG", "L");
		code.put("ATT", "I");
		code.put("ATC", "I");
		code.put("ATA", "I");
		code.put("ATG", "M");
		code.put("GTT", "V");
		code.put("GTC", "V");
		code.put("GTA", "V");
		code.put("GTG", "V");
		code.put("TCT", "S");
		code.put("TCC", "S");
		code.put("TCA", "S");
		code.put("TCG", "S");
		code.put("CCT", "P");
		code.put("CCC", "P");
		code.put("CCA", "P");
		code.put("CCG", "P");
		code.put("ACT", "T");
		code.put("ACC", "T");
		code.put("ACA", "T");
		code.put("ACG", "T");
		code.put("GCT", "A");
		code.put("GCC", "A");
		code.put("GCA", "A");
		code.put("GCG", "A");
		code.put("TAT", "Y");
		code.put("TAC", "Y");
		code.put("TAA", "*");
		code.put("TAG", "+");
		code.put("CAT", "H");
		code.put("CAC", "H");
		code.put("CAA", "Q");
		code.put("CAG", "Q");
		code.put("AAT", "N");
		code.put("AAC", "N");
		code.put("AAA", "K");
		code.put("AAG", "K");
		code.put("GAT", "D");
		code.put("GAC", "D");
		code.put("GAA", "E");
		code.put("GAG", "E");
		code.put("TGT", "C");
		code.put("TGC", "C");
		code.put("TGA", "#");
		code.put("TGG", "W");
		code.put("CGT", "R");
		code.put("CGC", "R");
		code.put("CGA", "R");
		code.put("CGG", "R");
		code.put("AGT", "S");
		code.put("AGC", "S");
		code.put("AGA", "R");
		code.put("AGG", "R");
		code.put("GGT", "G");
		code.put("GGC", "G");
		code.put("GGA", "G");
		code.put("GGG", "G");
		return code;
	}
	public  static ArrayList<BitSet> FastaToBitSet(ArrayList<Fasta> seqs){
		ArrayList<BitSet> newSeqs=new ArrayList<BitSet>();
		for(int i=0;i<seqs.size();i++){
			String sequence=seqs.get(i).getSequence().toUpperCase();
			if(sequence.contains("N")){
				System.err.println("Left out a "+sequence.length()+"bp region since it contained Ns.");
				continue;
			}
			newSeqs.add(DNAmanipulations.codeDNA(sequence));
		}
		return newSeqs;
	}
	public  static ArrayList<BitSet> toBitSet(ArrayList<String> seqs){
		ArrayList<BitSet> newSeqs=new ArrayList<BitSet>();
		for(int i=0;i<seqs.size();i++){
			if(seqs.get(i).toUpperCase().contains("N")){
				System.err.println("Left out a "+seqs.get(i).length()+"bp region since it contained Ns.");
				continue;
			}
			newSeqs.add(DNAmanipulations.codeDNA(seqs.get(i).toUpperCase()));
		}
		return newSeqs;
	}
	public static String decodeDNA(BitSet code){
		
		StringBuilder DNA=new StringBuilder("");
		for(int i=0;i<code.length()-1;i+=2){
			if(code.get(i)==false && code.get(i+1)==false){
				DNA.append('A');
			}else if(code.get(i)==false && code.get(i+1)==true){
				DNA.append('T');

			}else if(code.get(i)==true && code.get(i+1)==false){
				DNA.append('C');

			}else if(code.get(i)==true && code.get(i+1)==true){
				DNA.append('G');

			}
		}
		return DNA.toString();
	}
	

	public static String decodeDNA(long code,int length){
		StringBuilder DNA=new StringBuilder("");
		for(int i=0;i<length;i+=1){
			if(code%4==0){
				DNA.append('A');
			}else if(code%4==1){
				DNA.append('T');

			}else if(code%4==2){
				DNA.append('C');

			}else if(code%4==3){
				DNA.append('G');

			}
			code/=4;
		}
		return DNA.toString();
	}
	
	//for quick coding
	
	public static long codeDNALong(String dna){
		long code=0;
		long A=0;
		long T=1;
		long C=2;
		long G=3;
		for (int i=0;i<dna.length();i++){
			long newLetter=0;
			if(dna.charAt(i)=='A'){
				newLetter=A << i*2l;
				
			}else if(dna.charAt(i)=='T'){
				newLetter=T << i*2l;

			}else if(dna.charAt(i)=='C'){
				newLetter=C << i*2l;

			}else if(dna.charAt(i)=='G'){
				newLetter=G << i*2l;

			}else{
				System.err.println("Wrong letter in DNA sequence: "+dna.charAt(i));
				return -1;
			}
			code+=newLetter;
		}
		return code;
	}
	//for long coding...slow
//	public static BitSet codeDNA(String dna){
//		BitSet code=new BitSet(dna.length()*2);
//		for (int i=0;i<dna.length();i++){
//			if(dna.charAt(i)=='A'){
//				code.set(i*2,false);
//				code.set(i*2+1,false);
//			}else if(dna.charAt(i)=='T'){
//				code.set(i*2,false);
//				code.set(i*2+1,true);
//			}else if(dna.charAt(i)=='C'){
//				code.set(i*2,true);
//				code.set(i*2+1,false);
//			}else if(dna.charAt(i)=='G'){
//				code.set(i*2,true);
//				code.set(i*2+1,true);
//			}else{
//				System.err.println("Wrong letter in DNA sequence: "+dna.charAt(i));
//				return null;
//			}
//		}
//		return code;
//	}
	public static BitSet codeDNA(String dna){
		BitSet code=new BitSet(dna.length()*2);
		dna=dna.toUpperCase();
		for (int i=0;i<dna.length();i++){
			if(dna.charAt(i)=='A'){
				code.set(i*2,false);
				code.set(i*2+1,false);
			}else if(dna.charAt(i)=='T'){
				code.set(i*2,false);
				code.set(i*2+1,true);
			}else if(dna.charAt(i)=='C'){
				code.set(i*2,true);
				code.set(i*2+1,false);
			}else if(dna.charAt(i)=='G'){
				code.set(i*2,true);
				code.set(i*2+1,true);
			}else{
				System.err.println("Wrong letter in DNA sequence: "+dna.charAt(i));
				return null;
			}
		}
		code.set((dna.length())*2,true);
		return code;
	}
}
