package buildTreesFromSequenceData;

import java.io.*;
import java.util.*;

import util.phylogenetics.RunTreePrograms;




public class PerformBowtie {
	public static ArrayList<File> runMultiple(File core,ArrayList<File> cutSequences,File outFolder,String bowtiepath,String buildpath,boolean runBowtie,int seedLength,File bowtieOptions){
		ArrayList<File> alignFiles=new ArrayList<File>();
		for(int i=0;i<cutSequences.size();i++){
			
			
			String suffix=samtoolsExist()?".bam":".sam";
			File out=new File(outFolder+"/"+cutSequences.get(i).getName().split("\\.")[0]+suffix);
			runBowtie(core,cutSequences.get(i),out,bowtiepath,buildpath,runBowtie, seedLength,bowtieOptions);
			alignFiles.add(out);
		}
		return alignFiles;
	}
	
	
	public static boolean samtoolsExist(){
		String samtools="samtools";
		boolean bam=true;
		try{
			Runtime.getRuntime().exec(samtools);
			
		}catch(IOException e){
			//System.out.println("Samtools do not exist. Alignment output will be in SAM format.");
			bam=false;
		}
		if(bam){
			//System.out.println("Samtools do exist. Alignment output will be in BAM format.");
			
		}
		return bam;
	}
	
	public static void runBowtie(File core,File cut,File out,String bowtiepath,String buildpath,boolean runBowtie,int seedLength,File bowtieOptions){
		File database=new File(core+".1.bt2");
		//database creation
		try{
			if(!database.exists()){
				
				String buildDB=buildpath+" -f "+core+ " "+core;
				Process p=Runtime.getRuntime().exec(buildDB);
				if(p.waitFor()!=0){
					System.err.println(buildDB);
					InputStream i=p.getErrorStream();
					int c=0;
					while((c=i.read())!=-1){
						System.err.print((char)c);
					}
					System.err.println("Building the bowtie database was not successful!");
					System.exit(-1);
				}
			}
			//TODO
			HashMap<String,Boolean> paraHM=RunTreePrograms.getParameters(bowtieOptions);
			String paraLine=RunTreePrograms.getParametersLine(bowtieOptions);
			String seedString=paraHM.containsKey("-L")?"":" -L "+seedLength+" ";
			String N=paraHM.containsKey("-N")?"":" -N 1 ";
			String repeats=paraHM.containsKey("-k")?"":" -a ";
			//execute bowtie if it does not already exist
			if(!out.exists()||runBowtie){
				File temp=new File(out+".temp");
				
				String SAMCommand=bowtiepath+" -x "+core+" -U "+ cut+" -S "+out+ N+" --no-unal -a  --no-head "+seedString+paraLine;
				String BAMCommand=bowtiepath+" -x "+core+" -U "+ cut+" -S "+temp+ N+" --no-unal "+repeats+seedString+" "+paraLine;
				if(out.toString().endsWith(".sam")){
					
					runSAM(SAMCommand);
				}else{

					runBAM(BAMCommand,temp,out);
					temp.delete();
				}
				
			}
		}catch(InterruptedException e){
			e.printStackTrace();
			System.exit(-1);
		}catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
		
	}
	
	
	public static void runBAM(String bowtiecom,File temp,File out)throws InterruptedException,IOException{
		
		runSAM(bowtiecom);
		
		
		String samtoolsCom="samtools view -bS "+temp+" -o "+out;

		Process samtools=Runtime.getRuntime().exec(samtoolsCom);
		if(samtools.waitFor()!=0){
			System.err.println(bowtiecom);
			InputStream i=samtools.getErrorStream();
			int c=0;
			while((c=i.read())!=-1){
				System.err.print((char)c);
			}
			System.err.println("Bowtie was not successful!");
			System.exit(-1);
		}
		



	}
	
	public static void  runSAM(String bowtiecom)throws InterruptedException,IOException{
		Process p=Runtime.getRuntime().exec(bowtiecom);
		if(p.waitFor()!=0){
			System.err.println(bowtiecom);
			InputStream i=p.getErrorStream();
			int c=0;
			while((c=i.read())!=-1){
				System.err.print((char)c);
			}
			System.err.println("Bowtie was not successful!");
			System.exit(-1);
		}
	}
	
	public static void deleteDatabases(File folder){
		
		File[] list=folder.listFiles();
		for(int i=0;i<list.length;i++){
			if(list[i].getName().endsWith("bt2")){
				list[i].delete();
			}
		}

	}

}
