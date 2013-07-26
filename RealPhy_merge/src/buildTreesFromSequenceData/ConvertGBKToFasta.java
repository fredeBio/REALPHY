package buildTreesFromSequenceData;

import java.io.*;
import java.util.*;

import util.*;

public class ConvertGBKToFasta {
	public static void main(String args[]){
		File gbk=new File(args[0]);
		int flank=Integer.parseInt(args[2]);
		File out = new File(args[3]);
		writeCDS(gbk,out,flank,true);
	}
	
	public static HashMap<String,Fasta> writeCDS(File gbk,File out,int flank,boolean DNA){
		HashMap<String,Fasta> hm=new HashMap<String, Fasta>();
		try{
			if(!DNA){
				flank=0;
			}
			
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			ReadGenbank rgb=new ReadGenbank(gbk);
			HashMap<String,String> genome=Fasta.fasToHash(rgb.getSequence(),true);
			
			ArrayList<String> ids=rgb.getIds();
			int count=0;
			for(int j=0;j<ids.size();j++){
				String LocusId=ids.get(j);
				ArrayList<Info> CDS=rgb.getIntervals("CDS",LocusId);
				String seq=genome.get(LocusId);
				if(seq==null){
					System.err.println("Cannot find "+LocusId+" from genbank file "+gbk+ " in corresponding sequence file!");
					System.exit(-1);
				}
				//System.out.println(CDS.size());
				String contig=genome.size()>1?LocusId:"n/a";
				for(int i=0;i<CDS.size();i++){	
					String id=""+(count+1);
					
					String description=(count+1)+" "+CDS.get(i).getStart()+".."+CDS.get(i).getEnd()+" "+CDS.get(i).getOrient()+" "+contig+" "+CDS.get(i).getInfo();
					bw.write(">"+description+"\n");
					
					String subseq=getDNASequence(CDS.get(i),flank,seq);	
					if(DNA){
						hm.put(id, new Fasta(description,subseq));
						bw.write(subseq+"\n");
					}else{
						hm.put(id,new Fasta(description, DNAmanipulations.translate(subseq, DNAmanipulations.code())));
						bw.write(DNAmanipulations.translate(subseq, DNAmanipulations.code())+"\n");
					}
					count++;
				}
			}
			bw.close();
		}catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
		return hm;
	}
	public static HashMap<String,Fasta> writeCDS(File gbk,ArrayList<Fasta> genome,File out,int flank,boolean DNA){
		HashMap<String,Fasta> hm=new HashMap<String, Fasta>();
		try{
			if(!DNA){
				flank=0;
			}
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			ReadGenbank rgb=new ReadGenbank(gbk);
			int x=0;
			for(int j=0;j<genome.size();j++){
				ArrayList<Info> CDS=rgb.getIntervals("CDS",genome.get(j).getIdent().split("\\s+|>")[0]);
				String seq=genome.get(j).getSequence();
				for(int i=0;i<CDS.size();i++){	
					x++;
					String id=""+x;
					String description=x+" "+CDS.get(i).getStart()+".."+CDS.get(i).getEnd()+" "+CDS.get(i).getOrient()+" "+out.getName().split("\\.")[0]+" "+CDS.get(i).getInfo();
					bw.write(">"+description+"\n");
					String subseq=getDNASequence(CDS.get(i),flank,seq);	
					if(DNA){
						hm.put(id, new Fasta(description,subseq));
						bw.write(subseq+"\n");
					}else{
						hm.put(id,new Fasta(description, DNAmanipulations.translate(subseq, DNAmanipulations.code())));
						bw.write(DNAmanipulations.translate(subseq, DNAmanipulations.code())+"\n");
					}
				}
			}
			
			bw.close();
		}catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
		return hm;
	}
	public static ArrayList<Fasta> getCDS(File gbk,ArrayList<Fasta> genome,int flank,boolean DNA){
		ArrayList<Fasta> list=new ArrayList<Fasta>();
		if(!DNA){
			flank=0;
		}
		ReadGenbank rgb=new ReadGenbank(gbk);
		int x=0;
		for(int j=0;j<genome.size();j++){
			String seq=genome.get(j).getSequence();
			String id=genome.get(j).getIdent().split(">|\\s+")[0];
			ArrayList<Info> CDS=rgb.getIntervals("CDS",id);
			for(int i=0;i<CDS.size();i++){	
				x++;
				String description=x+" "+CDS.get(i).getStart()+".."+CDS.get(i).getEnd()+" "+CDS.get(i).getOrient()+" "+CDS.get(i).getInfo();
				String subseq=getDNASequence(CDS.get(i),flank,seq);	
				if(DNA){
					list.add(new Fasta(description,subseq));
				}else{
					list.add(new Fasta(description, DNAmanipulations.translate(subseq, DNAmanipulations.code())));
				}
			}
		}
		return list;
	}
	
	public static String getDNASequence(Info inf,int flank,String seq){
		boolean reverse=inf.info.endsWith("complement");
		//System.out.println(inf);
		int start=Math.min(inf.getStart()-1,inf.getEnd())-flank;
		int end=Math.max(inf.getStart()-1,inf.getEnd())+flank;
		String startAdd="";
		String endAdd="";
		if(start<0){
			startAdd=seq.substring(seq.length()+start);
			start=0;
		}
		if(end>seq.length()){
			endAdd=seq.substring(0,end-seq.length());
			end=seq.length()-1;
		}
		String subseq=seq.substring(start,end);
		subseq=startAdd+subseq+endAdd;
		subseq=reverse?DNAmanipulations.reverse(subseq):subseq;	
		return subseq;
	}
	
}
