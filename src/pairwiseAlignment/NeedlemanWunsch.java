package pairwiseAlignment;

import java.io.*;
import java.util.*;

import util.*;

public class NeedlemanWunsch {
	double dE,dO;
	double[][] A;
	double[][] B;
	double[][] C;
	Alignment al;
	HashMap<Character,HashMap<Character,Integer>> subMatrix=new HashMap<Character, HashMap<Character,Integer>>();
	String sA; String sB;


	
	
	
	public static HashMap<Character,HashMap<Character,Integer>> readSimilarityMatrix(File in){
		HashMap<Character,HashMap<Character,Integer>> hm=new HashMap<Character, HashMap<Character,Integer>>();

		try{
			BufferedReader br=new BufferedReader(new FileReader(in));
			String line="";
			boolean first=true;
			String[] order=null;
			int i=0;
			while((line=br.readLine())!=null){
				if(line.startsWith("#"))continue;
				
				String[] split=line.split("\\s+");
				if(first){
					first=false;
					order=split;
				}else{
					i++;
					char c=order[i].charAt(0);
					
					hm.put(c, new HashMap<Character, Integer>());
					for(int j=1;j<split.length;j++){
						hm.get(c).put(order[j].charAt(0), Integer.parseInt(split[j]));
					}
				}
				
			}
			br.close();
		}catch(IOException e){
			e.printStackTrace();
		}		return hm;
	}
	
	
	public static void main(String[] args)
	{
		//	        int[] ar = NeedlemanWunsch.convertStringToArr("TCGAAGCT");        
		//	        System.out.println(ar);
		File input=new File(args[0]);
		File in=new File(args[1]);
		ArrayList<Fasta> fas=Fasta.readFasta(input);
		String a=fas.get(0).getSequence().toUpperCase();
		String b=fas.get(1).getSequence().toUpperCase();
		NeedlemanWunsch nw=new NeedlemanWunsch(a, b,readSimilarityMatrix(in),10,.1);
//		for (int y = 0; y < ar.length; y++)
//		{
//			System.out.println("");
//			for (int x = 0; x < ar[y].length; x++)
//				System.out.print(ar[y][x] +"\t ");
//		}
//		System.out.println();
		
		System.out.println(nw.getPairwiseIdentity());
		System.out.println(nw.getAlignments());
	}
	public NeedlemanWunsch(String a,String b,HashMap<Character,HashMap<Character,Integer>> matrix,double gapOpen,double gapExtend){
		a=checkForStops(a);
		b=checkForStops(b);
		sA=a;
		sB=b;
		subMatrix=matrix;
		dE=gapExtend;
		dO=gapOpen;
		
		calculateMatrix();
		doAlignments();
	}
	public String checkForStops(String a){
		StringBuffer temp=new StringBuffer();
		for(int i=0;i<a.length();i++){
			if(a.charAt(i)=='#'||a.charAt(i)=='*'||a.charAt(i)=='+'){
				
			}else{
				temp.append(a.charAt(i));
			}
		}
		return temp.toString();
	}
	
	public Alignment getAlignments(){
		return al;
	}
	
	public  double getScore(){
		return max(A[A.length-1][A[A.length-1].length-1],B[B.length-1][B[B.length-1].length-1],C[C.length-1][C[C.length-1].length-1]);
	}
	
	public  static double getPairwiseIdentity(String a,String b){
		if(a.length()!=b.length()){
			System.err.println("String a has not the same length as string b. Return NaN.");
			return Double.NaN;
		}
		int identities=0;
		for(int i=0;i<a.length();i++){
			if(a.charAt(i)==b.charAt(i)){
				identities++;
			}
		}
		return (1.0*identities)/((a.length()+b.length())/2);
		
	}
	
	public  static double getPairwiseDifferences(String a,String b){
		if(a.length()!=b.length()){
			System.err.println("String a has not the same length as string b. Return NaN.");
			return Double.NaN;
		}
		int differences=0;
		for(int i=0;i<a.length();i++){
			if(a.charAt(i)!=b.charAt(i)){
				differences++;
			}
		}
		return (1.0*differences)/((a.length()+b.length())/2);
		
	}
	
	public  double getPairwiseIdentity(){
		int identities=0;
		for(int i=0;i<al.alA.length();i++){
			if(al.alA.charAt(i)==al.alB.charAt(i)){
				identities++;
			}
		}
		return (1.0*identities)/((al.alA.length()+al.alB.length())/2);
		
	}
	
	public  void doAlignments()
	{
		sA=sA.toUpperCase();
		sB=sB.toUpperCase();
		//int[] A=convertStringToArr(sA); int[] B=convertStringToArr(sB);
		String alA = "";
		String alB = "";        
		int i = sA.length();
		int j = sB.length();
//		for (int y = 0; y < M.length; y++)
//		{
//			System.out.println("");
//			for (int x = 0; x < M[y].length; x++)
//				System.out.print(M[y][x] +"\t ");
//		}
//		System.out.println();
//		
//		for (int y = 0; y < I.length; y++)
//		{
//			System.out.println("");
//			for (int x = 0; x < I[y].length; x++)
//				System.out.print(I[y][x] +"\t ");
//		}
//		System.out.println();
		double score=max(A[i][j],B[i][j],C[i][j]);
		char mat=A[i][j]==score?'a':B[i][j]==score?'b':'c';
		while (i > 0 && j > 0)
		{						
			double scoreleft=B[i-1][j];
			double scoreup=C[i][j-1];
			if (mat=='a')
			{
				score=max(A[i-1][j-1],B[i-1][j-1],C[i-1][j-1]);
				mat=A[i-1][j-1]==score?'a':B[i-1][j-1]==score?'b':'c';
				alA = sA.charAt(i-1) + alA;
				alB = sB.charAt(j-1) + alB;
				i--;j--;                
			}
			else if (mat=='b')
			{
				score=max(A[i-1][j]-dO,B[i-1][j]-dE,C[i-1][j]-dO);
				mat=A[i-1][j]-dO==score?'a':B[i-1][j]-dE==score?'b':'c';
				alA = sA.charAt(i-1) + alA;
				alB = "-" + alB;
				i--;
			}
			else if(mat=='c')
			{
				score=max(A[i][j-1]-dO,B[i][j-1]-dO,C[i][j-1]-dE);
				mat=A[i][j-1]-dO==score?'a':B[i][j-1]-dO==score?'b':'c';
				
				alA = "-" + alA;
				alB = sB.charAt(j-1) + alB;
				j--;
			}else{
				System.err.println("Error: "+score+" "+scoreup+" "+scoreleft+" ");
				break;
			}
		}
		while(i > 0)
		{
			alA = sA.charAt(i - 1) + alA;
			alB = "-" + alB;
			i--;            
		}
		while(j > 0)
		{
			alA = "-" + alA;
			alB = sB.charAt(j - 1) + alB;
			j--;            
		}
		al= new Alignment(alA,alB);
	}

	
	
	public  void calculateMatrix()
	{
		String source=sA.toUpperCase();
		String dest=sB.toUpperCase();
		A = new double[source.length()+1][dest.length()+1];
		B = new double[source.length()+1][dest.length()+1];
		C = new double[source.length()+1][dest.length()+1];

		for (int y = 0; y < source.length(); y++){
			A[y][0] = 0;
			B[y][0] = 0;
			C[y][0] = 0;
		}
		for (int x = 0; x < dest.length(); x++){
			A[0][x] =0;
			B[0][x]=0;
			C[0][x]=0;
		}
		for (int y = 1; y < source.length() + 1; y++){
			for (int x = 1; x < dest.length() +1; x++)
			{                    
				double k = A[y-1][x-1] + similar(source.charAt(y-1) , dest.charAt(x-1));
				double l = B[y-1][x-1] + similar(source.charAt(y-1) , dest.charAt(x-1));
				double m = C[y-1][x-1] + similar(source.charAt(y-1) , dest.charAt(x-1));
				A[y][x]=max(k,l,m);
				k=A[y-1][x]-dO;
				l=B[y-1][x]-dE;
				m=C[y-1][x]-dO;
				B[y][x]=max(k,l,m);
				k=A[y][x-1]-dO;
				l=B[y][x-1]-dO;
				m=C[y][x-1]-dE;
				C[y][x] = max(k,l,m);

			}

		}
		
	}
	double max(double... paras){
		double max=Double.NEGATIVE_INFINITY;
		for(int i=0;i<paras.length;i++){
			if(paras[i]>max){
				max=paras[i];
			}
		}
		return max;
	}
	public  int similar(char first, char second)
	{
		//System.out.println(first+" "+second);
		return subMatrix.get(first).get(second);
	}


}


