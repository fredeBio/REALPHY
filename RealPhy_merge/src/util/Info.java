package util;

import java.io.*;
import java.util.*;


public class Info extends Interval<Info> {
		/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

		public Info(int S,int E,String inf,int RepNumber){
			super(S,E);
			info=inf;
			repNumber=RepNumber;
		}
//		public Info(int S,int E,String inf,char Orientation){
//			super(S,E);
//			info=inf;
//			repNumber=-1;
//			orient=Orientation;
//		}
//		public Info(int S,int E,String inf,char Orientation,int Codon_start){
//			super(S,E);
//			info=inf;
//			repNumber=-1;
//			orient=Orientation;
//			codon_start=Codon_start;
//		}
		public Info(int S,int E,String inf){
			super(S,E);
			info=inf;
			repNumber=-1;
			orient='?';
		}
		public Info(Info inf){
			super(inf.start,inf.end);
			info=inf.info;
			orient=inf.orient;
			pseudo=inf.pseudo;
			codon_start=inf.codon_start;
			repNumber=inf.repNumber;
			feature=inf.feature;
					

		}
		public String info;
		public int repNumber=-1;
		char orient;
		int codon_start;
		boolean pseudo;
		String feature;
		
		public Info setFeature(String Feature){
			feature=Feature;
			return this;
		}
		
		public String getFeature(){
			return feature;
		}
		
		public char getOrient(){
			return orient;
		}
		public int getCodonStart(){
			return codon_start;
		}
		public Info setCodonStart(int cs){
			codon_start=cs;
			return this;
		}
		
		public boolean getPseudo(){
			return pseudo;
		}
		
		public Info setPseudo(boolean Pseudo){
			pseudo=Pseudo;
			return this;
		}
		
		public void append(Info s){
			info+=" "+s.info;
		}

		public String getInfo(){
			return info;
		}
		public void setStart(int s){
			start=s;
		}
		public int size(){
			return end-start+1;
		}
		public void setEnd(int e){
			end=e;
		}
		public String toString(){
			return start+"\t"+end+"\t"+info+"\t"+orient;
		}
		
		public String getGeneOrLocusInfo(){
			if(info.contains("gene= ")){
				return info.split("gene= ")[1].split("\\s+")[0];
			}else if(info.contains("locus_tag= ")){
				return info.split("locus_tag= ")[1].split("\\s+")[0];
			}
			
			return "";
		}
		public Info setOrient(char o){
			orient=o;
			return this;
		}
		public static void write(ArrayList<Info> infos,File out){
			try{
				BufferedWriter bw=new BufferedWriter(new FileWriter(out));
				for(int i=0;i<infos.size();i++){
					bw.write(infos.get(i)+"\n");
				}
				bw.close();
			}catch(IOException e){
				e.printStackTrace();
			}
		}
		
}
