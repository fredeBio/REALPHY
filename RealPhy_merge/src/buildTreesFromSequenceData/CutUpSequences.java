package buildTreesFromSequenceData;

import java.io.*;
import java.util.*;

import util.*;

public class CutUpSequences {
	public static void main(String args[]){
		runCutUpSequences(args);
	}
	
	public static ArrayList<File> runCutUpSequences(String args[]){
		ArrayList<File> cutSeqs=new ArrayList<File>();
		File inputFolder=new File(args[0]);
		int length=Integer.parseInt(args[1]);
		File outputFolder=new File(args[2]);
		boolean clean=new Boolean(args[3]);
		File[] files=inputFolder.listFiles();
		HashMap<String,Boolean> hm=new HashMap<String,Boolean>();
		for(int i=0;i<files.length;i++){
			String name=files[i].getName();
			
			if(RealPhy.hasExtension(name, RealPhy.fasExt)||RealPhy.hasExtension(name, RealPhy.gbkExt)){
				File temp=cutFasta(files[i],length,outputFolder,clean);
				
				if(temp!=null&&!hm.containsKey(temp.getName())){
					hm.put(temp.getName(),true);
					cutSeqs.add(temp);
				}
			}else if(files[i].getName().endsWith("fastq")){
				//CUT FASTQ FILES!!!
				
				
				File temp=cutFastq(files[i],length,outputFolder,clean);
				if(temp!=null&&!hm.containsKey(temp.getName())){
					hm.put(temp.getName(),true);
					cutSeqs.add(temp);
				}
			}else if(files[i].getName().endsWith("fastq.gz")){
				cutSeqs.add(files[i]);

			}
		}
		return cutSeqs;
	}
	
	public static File cutFastaSimple(File fas,int length,File cutFolder,boolean clean){
		try{
			File idF=new File(fas.getName());
			String id=idF.getName().split("\\.")[0];
			if(!cutFolder.exists()){
				cutFolder.mkdir();
			}
			File out=new File(cutFolder+"/"+id+"_"+length+".fastq");
			if(out.exists()&&!clean){
				System.out.println(out+" already exists. Continue with next file.");
				return out;
			}
			BufferedReader br=new BufferedReader(new FileReader(fas));
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			String line="";
			String[] lines=new String[4];
			int k=0;
			Fasta fasta=new Fasta("","");
			int warning=0;
			while((line=br.readLine())!=null){
				
				if(k>0&&line.startsWith(">")){
					lines[0]="@"+fasta.getIdent();
					lines[1]=fasta.getSequence();
					lines[2]="+";
					StringBuffer sb=new StringBuffer();
					for(int j=0;j<lines[1].length();j++){
						sb.append('J');
					}
					lines[3]=sb.toString();
					if (lines[1].length()<length){
						warning++;
						if(warning==1)System.err.println("In file "+fas.getName()+" at least one of the fasta sequences is shorter ("+lines[1].length()+") than the specified readLength("+length+")!");
					}
					write(lines,length,bw);
					fasta=new Fasta(line.trim().substring(1),"");
				}else if(!line.startsWith(">")){
					line=line.replaceAll("\\s+", "");
					fasta.setSequence(fasta.getSequence()+line);
					
				}else if(k==0&&line.startsWith(">")){
					fasta.setIdent(line.trim().substring(1));
				}
				k=1;
			}
			lines[0]="@"+fasta.getIdent();
			lines[1]=fasta.getSequence();
			lines[2]="+";
			StringBuffer sb=new StringBuffer();
			for(int j=0;j<lines[1].length();j++){
				sb.append('J');
			}
			lines[3]=sb.toString();
			write(lines,length,bw);
			bw.close();
			br.close();
			return out;
		}catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
			return null;
		}
	}

	
	public static File cutFastq(File fastq,int length,File cutFolder,boolean clean){
		try{
			File idF=new File(fastq.getName());
			String id=idF.getName().split("\\.")[0];
			if(!cutFolder.exists()){
				cutFolder.mkdir();
			}
			File out=new File(cutFolder+"/"+id+"_"+length+".fastq");
			if(out.exists()&&!clean){
				System.out.println(out+" already exists. Continue with next file.");
				return out;
			}
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			BufferedReader br=new BufferedReader(new FileReader(fastq));
			String line="";
			int j=0;
			String[] lines=new String[4];
			int warning=0;
			boolean write=false;
			while((line=br.readLine())!=null){
				line=line.trim();
				if(line.length()==0){
					continue;
				}
				if(j%4==0&&j>0){
					if (lines[1].length()<length){
						warning++;
						if(warning==1)System.err.println("In file "+fastq.getName()+" at least one of the fastq sequences is shorter ("+lines[1].length()+") than the specified readLength("+length+")!");
					}
					if(lines[1].length()>30){
						write(lines,length,bw);
						write=true;
					}
				}
				lines[j%4]=line;	
				j++;
			}
			write(lines,length,bw);
			br.close();
			bw.close();
			if(write){
				return out;
			}else return null;
		}catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
			return null;
		}
	}
	
	public static void write(String[] lines,int splitLength,BufferedWriter bw) throws IOException {
		String id=lines[0];
		String seq=lines[1];
		String quality=lines[3];
		if (seq.length()<splitLength){
			bw.write(id+"_0\n"+seq+"\n"+lines[2]+"\n"+quality+"\n");

		}else{
			for(int i=0;i<seq.length()/splitLength;i++){
				bw.write(id+"_"+i+"\n"+seq.substring(i*splitLength,(i+1)*splitLength)+"\n"+lines[2]+"\n"+quality.substring(i*splitLength,(i+1)*splitLength)+"\n");
			}
		}
		
		
	}

	/**
	 * Cuts up a fasta file into "length" pieces and stores them in a fastq file. Each base quality is set to the J PHRED score. 
	 * The reads are named after the nucleotide position in the FASTA file.
	 * @param fasta
	 * Fasta file that will be cut.
	 * @param length
	 * Length of short reads that are being produced.
	 * @param cutFolder
	 * The folder in which the newly created read file will be stored.
	 * @param clean
	 * If true, then old cut files in the cut folder will be overwritten.
	 * @return
	 * Returns the location of the newly created fastq file.
	 */
	
	public static File cutFasta(File fasta,int length,File cutFolder,boolean clean){
		try{
			File idF=new File(fasta.getName());
			String id=idF.getName().split("\\.")[0];
			if(!cutFolder.exists()){
				cutFolder.mkdir();
			}
			File out=new File(cutFolder+"/"+id+"_"+length+"fasta.fastq");
			if(out.exists()&&!clean){
				System.out.println(out+" already exists. Continue with next file.");
				return out;
			}
			String name=fasta.getName();
			ArrayList<Fasta> fas=RealPhy.hasExtension(name, RealPhy.fasExt)?Fasta.readFasta(fasta):new ReadGenbank(fasta).getSequence();
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			boolean write=false;
			int size=fas.size();
			//counts the total number of nucleotides in the fasta file
			//for merging purposes, hence the name of the individual contigs will be irrelevant
			int total=0;
			for(int j=0;j<size;j++){
				String seq=fas.get(j).getSequence();
				int seqLength=seq.length();
				for(int k=0;k<=seqLength-length;k++){
					write=true;
					bw.write("@"+total+"\n");
					bw.write(seq.substring(k,k+length)+"\n");
					bw.write("+\n");
					total++;
					for(int m=0;m<length;m++)
						bw.write("J");
					bw.write("\n");
				}
			}
			
			bw.close();
			if(!write){
				out.delete();
				return null;
			}else return out;
		}catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
		return null;
	}
	
}
