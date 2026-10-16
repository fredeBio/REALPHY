package util;

import java.io.*;
import java.util.*;

public class FastQ {
	String line1="";
	String line2="";
	String line3="";
	String line4="";
	public FastQ(String l1,String l2,String l3,String l4){
		line1=l1;
		line2=l2;
		line3=l3;
		line4=l4;
	}
	public FastQ(){
	}
	public void setLine1(String l1){
		line1=l1;
	}
	public void setLine2(String l2){
		line2=l2;
	}
	public void setLine3(String l3){
		line3=l3;
	}
	public void setLine4(String l4){
		line4=l4;
	}

	public static ArrayList<Fasta> toFasta(ArrayList<FastQ> fastq){
		ArrayList<Fasta> fasta=new ArrayList<Fasta>();
		for(int i=0;i<fastq.size();i++){
			fasta.add(new Fasta(fastq.get(i).line1,fastq.get(i).line2));
		}
		return fasta;
	}
	
	public static int getQualityPos(String line,char q){
		for(int i=0;i<line.length();i++){
			if(line.charAt(i)==q){
				return i;
			}
		}
		return line.length();
	}	
	
	public static PairedEnd FilterPairedEnd(PairedEnd pe,char quality,int minSeqSize,String identAddition){
		PairedEnd filtered=new PairedEnd(new ArrayList<FastQ>(),new ArrayList<FastQ>());
		int avgLength=0;
		for(int i=0;i<pe.p1.size();i++){
			int q1=getQualityPos(pe.p1.get(i).line4,quality);
			int q2=getQualityPos(pe.p2.get(i).line4,quality);
			int q=q1>q2?q2:q1;
			if(q>minSeqSize){
				filtered.p1.add(pe.p1.get(i));
				filtered.p2.add(pe.p2.get(i));
				avgLength+=q;
			}
		}
		System.out.println(avgLength/(filtered.p1.size()*1.0));
		return filtered;
	}
	public static void write(ArrayList<FastQ> fastq,File out,boolean append){
		try{
			
			BufferedWriter bw=new BufferedWriter(new FileWriter(out,append));
			for(int i=0;i<fastq.size();i++){
				bw.write(fastq.get(i).line1+"\n"+fastq.get(i).line2+"\n"+fastq.get(i).line3+"\n"+fastq.get(i).line4+"\n");
			}
			bw.close();
		}catch(IOException e){
			e.printStackTrace();
		}
	}
	public static void writePairedEnd(PairedEnd pe,File out1,File out2, boolean append){
		write(pe.p1,out1,append);
		write(pe.p2,out2,append);
	}
	
	public static ArrayList<Fasta> toFasta(ArrayList<FastQ> fastq,char quality,int minSeqSize,String identAddition){
		ArrayList<Fasta> fasta=new ArrayList<Fasta>();
		int avgLength=0;
		for(int i=0;i<fastq.size();i++){
			int q=getQualityPos(fastq.get(i).line4,quality);
			if(q>minSeqSize){
				fasta.add(new Fasta(fastq.get(i).line1+identAddition,fastq.get(i).line2.substring(0,q)));
				avgLength+=q;
			}else{
				fasta.add(new Fasta(fastq.get(i).line1+identAddition,"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"));
			}
		}
		System.out.println(avgLength/(fasta.size()*1.0));
		return fasta;
	}
	public static ArrayList<Fasta> toFasta(ArrayList<FastQ> fastq,char quality){
		ArrayList<Fasta> fasta=new ArrayList<Fasta>();
		int avgLength=0;
		for(int i=0;i<fastq.size();i++){
			int q=getQualityPos(fastq.get(i).line4,quality);
			fasta.add(new Fasta(fastq.get(i).line1,fastq.get(i).line2.substring(0,q)));
			avgLength+=q;
		}
		System.out.println(avgLength/(fasta.size()*1.0));
		return fasta;
	}
	
	public static void ReadAndFilterPairedEnd(File in1,File in2,File out1,File out2, char quality,int minSeqSize,int buffersize,String identAddition){
		PairedEnd pe=new PairedEnd(new ArrayList<FastQ>(),new ArrayList<FastQ>());
		try{
			out1.delete();
			out2.delete();
			BufferedReader br1=new BufferedReader(new FileReader(in1));
			BufferedReader br2=new BufferedReader(new FileReader(in2));
			String line1="";
			String line2="";
			int lineNumber=0;
			while((line1=br1.readLine())!=null && (line2=br2.readLine())!=null){
				if(lineNumber%4==0){
					if((lineNumber/4)%buffersize==0 && lineNumber>0){
						PairedEnd filter=FilterPairedEnd(pe, quality, minSeqSize, identAddition);
						writePairedEnd(filter, out1, out2,true);
						pe=new PairedEnd(new ArrayList<FastQ>(),new ArrayList<FastQ>());
					}
					FastQ newEntry1=new FastQ();
					newEntry1.setLine1(line1);
					pe.p1.add(newEntry1);
					
					FastQ newEntry2=new FastQ();
					newEntry2.setLine1(line2);
					pe.p2.add(newEntry2);
				}else if(lineNumber%4==1){
					pe.p1.get(pe.p1.size()-1).setLine2(line1);
					pe.p2.get(pe.p2.size()-1).setLine2(line2);
				}else if(lineNumber%4==2){
					pe.p1.get(pe.p1.size()-1).setLine3(line1);
					pe.p2.get(pe.p2.size()-1).setLine3(line2);
				}else if(lineNumber%4==3){
					pe.p1.get(pe.p1.size()-1).setLine4(line1);
					pe.p2.get(pe.p2.size()-1).setLine4(line2);
				}
				lineNumber++;
			}
			br1.close();
			br2.close();
			PairedEnd filter=FilterPairedEnd(pe, quality, minSeqSize, identAddition);
			writePairedEnd(filter, out1, out2,true);
		}catch(IOException e){
			e.printStackTrace();
		}
	}
	
	public static void ReadFastQ(File in,File out, char quality,int minSeqSize,int buffersize,String identAddition){
		ArrayList<FastQ> al=new ArrayList<FastQ>();
		try{
			out.delete();
			BufferedReader br=new BufferedReader(new FileReader(in));
			String line="";
			int lineNumber=0;
			while((line=br.readLine())!=null){
				if(lineNumber%4==0){
					if((lineNumber/4)%buffersize==0 && lineNumber>0){
						ArrayList<Fasta> fas=toFasta(al, quality,minSeqSize,identAddition);
						Fasta.write(fas, out,true);
						al=new ArrayList<FastQ>();
					}
					
					
					FastQ newEntry=new FastQ();
					newEntry.setLine1(line);
					al.add(newEntry);
				}else if(lineNumber%4==1){
					al.get(al.size()-1).setLine2(line);
				}else if(lineNumber%4==2){
					al.get(al.size()-1).setLine3(line);
				}else if(lineNumber%4==3){
					al.get(al.size()-1).setLine4(line);
				}
				lineNumber++;
			}
			br.close();
			ArrayList<Fasta> fas=toFasta(al, quality,minSeqSize,identAddition);
			Fasta.write(fas, out,true);
			al=new ArrayList<FastQ>();
		}catch(IOException e){
			e.printStackTrace();
		}
	}
	public static ArrayList<FastQ> ReadFastQ(File in){
		ArrayList<FastQ> al=new ArrayList<FastQ>();
		try{
			BufferedReader br=new BufferedReader(new FileReader(in));
			String line="";
			int lineNumber=0;
			while((line=br.readLine())!=null){
				if(lineNumber%4==0){
					FastQ newEntry=new FastQ();
					newEntry.setLine1(line);
					al.add(newEntry);
				}else if(lineNumber%4==1){
					al.get(al.size()-1).setLine2(line);
				}else if(lineNumber%4==2){
					al.get(al.size()-1).setLine3(line);
				}else if(lineNumber%4==3){
					al.get(al.size()-1).setLine4(line);
				}
				lineNumber++;
			}


			br.close();
		}catch(IOException e){
			e.printStackTrace();
		}
		return al;
	}
	public static class PairedEnd{
		ArrayList<FastQ> p1;
		ArrayList<FastQ> p2;
		public PairedEnd(ArrayList<FastQ> P1,ArrayList<FastQ> P2){
			p1=P1;
			p2=P2;
		}
	}
	
}
