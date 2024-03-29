package solexa;

import java.io.Serializable;
import java.util.*;


public class Arrays implements Serializable{

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	double[][] array;
	double[] cov;
	HashMap<Integer,Integer>[] queryID;	
	public final int length;
	String ref;
	public Arrays(int Length){
		length=Length;
		array=new double [length][10];
		queryID=new HashMap[length];
	}
	public Arrays(String ref){
		length=ref.length()+1;
		this.ref=ref;
		array=new double [length][10];
		queryID=new HashMap[length];
	}
	public void setCoverage(double[] Coverage){
		cov=Coverage;
	}

	public double[] getCoverage(){
		return cov;
	}
//	public void set(int pos,String base){
//
//		if(base.equalsIgnoreCase("AF") && pos<array.length)array[pos][0]++;
//		else if(base.equalsIgnoreCase("TF") && pos<array.length)array[pos][1]++;
//		else if(base.equalsIgnoreCase("CF") && pos<array.length)array[pos][2]++;
//		else if(base.equalsIgnoreCase("GF") && pos<array.length)array[pos][3]++;
//		else if(base.equalsIgnoreCase("AR") && pos<array.length)array[pos][4]++;
//		else if(base.equalsIgnoreCase("TR") && pos<array.length)array[pos][5]++;
//		else if(base.equalsIgnoreCase("CR") && pos<array.length)array[pos][6]++;
//		else if(base.equalsIgnoreCase("GR") && pos<array.length)array[pos][7]++;
//	}
	
	
	public void set(int pos,String base,double weight){
		
		if(base.charAt(1)=='F'){
			if(base.charAt(0)=='A' && pos<array.length){array[pos][0]+=weight;return;}
			else if(base.charAt(0)=='T' && pos<array.length){array[pos][1]+=weight;return;}
			else if(base.charAt(0)=='C'&& pos<array.length){array[pos][2]+=weight;return;}
			else if(base.charAt(0)=='G' && pos<array.length){array[pos][3]+=weight;return;}
			else if(base.charAt(0)=='-' && pos<array.length){array[pos][4]+=weight;return;}

		}else if(base.charAt(1)=='R'){
			if(base.charAt(0)=='A' && pos<array.length){array[pos][5]+=weight;return;}
			else if(base.charAt(0)=='T' && pos<array.length){array[pos][6]+=weight;return;}
			else if(base.charAt(0)=='C'&& pos<array.length){array[pos][7]+=weight;return;}
			else if(base.charAt(0)=='G' && pos<array.length){array[pos][8]+=weight;return;}
			else if(base.charAt(0)=='-' && pos<array.length){array[pos][9]+=weight;return;}

		}
		//System.err.println("Do not recognize "+ base+"!");
		
	}
	
	public void set(int pos,String base,double weight,String readSequence){
		base=base.toUpperCase();
		if(ref!=null){
			char refb=ref.charAt(pos-1);
			if(complement.containsKey(refb)){
				char refc=complement.get(refb);
				if((refb+"F").equals(base)||(refc+"R").equals(base)){
					throw new RuntimeException("Substitutions that are identical with the reference are not possible!\n"+readSequence+"\n"+pos+"\n"+base+"\n"+ref);
				}
			}
		}
		set(pos,base,weight);
	}
	public static String[] bases=new String[]{"AF","TF","CF","GF","-F","AR","TR","CR","GR","-R"};

	static HashMap<Character,Character> complement=new HashMap<Character, Character>();
	static{
		complement.put('A', 'T');
		complement.put('T', 'A');
		complement.put('C', 'G');
		complement.put('G', 'C');
	}
	
	public void addPseudoCount(double count,String ref){
		for(int i=1;i<array.length;i++){
			cov[i]+=count*8;
			char r=ref.charAt(i-1);
			char c=complement.get(r);
			for(int j=0;j<array[i].length;j++){
				if((bases[j].charAt(1)=='F'&&r!=bases[j].charAt(0))||(bases[j].charAt(1)=='R'&&c!=bases[j].charAt(0))){
					array[i][j]+=count;
				}
			}
		}
	}
	public void set(int pos,String base,int readPos,double weight,String readSeq){
		set(pos,base,weight,readSeq);
		if(queryID[pos]==null){
			HashMap<Integer,Integer> temp=new HashMap<Integer,Integer>();
			temp.put(readPos,1);
			queryID[pos]=temp;
		}else{
			if(queryID[pos].containsKey(readPos)){
				int item=queryID[pos].get(readPos);
				item+=1;
				queryID[pos].put(readPos,item);
			}
		}
	}
	
	public boolean isset(int pos){
		if(queryID[pos]!=null){
			return true;
		}else return false;
		
	}
	
	public double numBases(int pos){
		double sum=0;
		for(int i=0;i<array[pos].length;i++){
			sum+=array[pos][i];
		}
		return sum;
	}
	
	public double numBasesNoGaps(int pos){
		double sum=0;
		for(int i=0;i<array[pos].length;i++){
			sum+=i!=9&&i!=4?array[pos][i]:0;
		}
		return sum;
	}
	
	public double get(int pos,String base){
		base=base.toUpperCase();
		if(base.charAt(1)=='F'){
			if(base.charAt(0)=='A' && pos<array.length)return array[pos][0];
			else if(base.charAt(0)=='T' && pos<array.length)return array[pos][1];
			else if(base.charAt(0)=='C'&& pos<array.length)return array[pos][2];
			else if(base.charAt(0)=='G' && pos<array.length)return array[pos][3];
			else if(base.charAt(0)=='-' && pos<array.length)return array[pos][4];

		}else if(base.charAt(1)=='R'){
			if(base.charAt(0)=='A' && pos<array.length)return array[pos][5];
			else if(base.charAt(0)=='T' && pos<array.length)return array[pos][6];
			else if(base.charAt(0)=='C'&& pos<array.length)return array[pos][7];
			else if(base.charAt(0)=='G' && pos<array.length)return array[pos][8];
			else if(base.charAt(0)=='-' && pos<array.length)return array[pos][9];

		}

		//System.err.println("Do not recognize "+ base+"!");
		return -1;

	}
	public double getFrequency(int pos){
		return (numBases(pos)*1.0)/(cov[pos]*1.0);
	}
	
	public double getFrequencyNoGaps(int pos){
		return (numBasesNoGaps(pos)*1.0)/(cov[pos]*1.0);
	}
	
	public HashMap<Integer,Integer> getQueryID(int pos){
		return queryID[pos];
	}
	public int getMaxNuc(int pos){
		double max=-1;
		int nuc=-1;
		for(int i=0;i<array[pos].length;i++){
			nuc=max>array[pos][i]?nuc:i;
			max=max>array[pos][i]?max:array[pos][i];
		}
		return nuc;
		
	}


	
	public char getMaxNucForward(int pos){
		int nuc=getMaxNuc(pos);
		
		return nuc==0||nuc==6?'A':nuc==1||nuc==5?'T':nuc==2||nuc==8?'C':nuc==3||nuc==7?'G':'-';
		
	}
	
	public double numPolymorphism(int pos,char base){
		char nuc=base;
		if(nuc=='A'){
			return array[pos][0]+array[pos][6];
		}else if(nuc=='T'){
			return array[pos][1]+array[pos][5];
		}else if(nuc=='C'){
			return array[pos][2]+array[pos][8];
		}else if(nuc=='G'){
			return array[pos][3]+array[pos][7];
		}else {
			return array[pos][4]+array[pos][9];
		}
	}
	
	public double numMajorPolymorphism(int pos){
		char nuc=getMaxNucForward(pos);
		if(nuc=='A'){
			return array[pos][0]+array[pos][6];
		}else if(nuc=='T'){
			return array[pos][1]+array[pos][5];
		}else if(nuc=='C'){
			return array[pos][2]+array[pos][8];
		}else if(nuc=='G'){
			return array[pos][3]+array[pos][7];
		}else {
			return array[pos][4]+array[pos][9];
		}
	}
	
//	private int getMax(int pos){
//		if(arrayA[pos]>=arrayT[pos] && arrayA[pos]>=arrayC[pos]&& arrayA[pos]>=arrayG[pos])return arrayA[pos];
//		else if(arrayT[pos]>=arrayA[pos] && arrayT[pos]>=arrayG[pos]&& arrayT[pos]>=arrayC[pos])return arrayT[pos];
//		else if(arrayC[pos]>=arrayA[pos] && arrayC[pos]>=arrayG[pos]&& arrayC[pos]>=arrayT[pos])return arrayC[pos];
//		return arrayG[pos];
//	}
	
}
