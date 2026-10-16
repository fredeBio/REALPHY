package buildTreesFromSequenceData;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;


public class PerformSoap {
	public static ArrayList<File> runMultipleSoaps(File core,ArrayList<File> cutSequences,File outFolder,String soappath,String buildpath,boolean runSoap){
		ArrayList<File> soapFiles=new ArrayList<File>();
		for(int i=0;i<cutSequences.size();i++){
			File out=new File(outFolder+"/"+cutSequences.get(i).getName().split("\\.")[0]+".sop");
			runSoap(core,cutSequences.get(i),out,soappath,buildpath,runSoap);
			soapFiles.add(out);
		}
		return soapFiles;
	}
	
	public static void runSoap(File core,File cut,File out,String soappath,String buildpath,boolean runSoap){
		File database=new File(core+".index.amb");
		//database creation
		try{
			if(!database.exists()){
				String buildDB=buildpath+" "+core;
				Process p=Runtime.getRuntime().exec(buildDB);
				if(p.waitFor()!=0){
					System.err.println(buildDB);
					InputStream i=p.getErrorStream();
					int c=0;
					while((c=i.read())!=-1){
						System.err.print((char)c);
					}
					System.err.println("Building the soap database was not successful!");
					System.exit(-1);
				}
			}
			//execute soap if it does not already exist
			if(!out.exists()||runSoap){
				String soapcom=soappath+" -o "+out+" -D "+core+".index -a "+cut+" -r 0";
				Process p=Runtime.getRuntime().exec(soapcom);
				if(p.waitFor()!=0){
					System.err.println(soapcom);
					InputStream i=p.getErrorStream();
					int c=0;
					while((c=i.read())!=-1){
						System.err.print((char)c);
					}
					System.err.println("Soap was not successful!");
					System.exit(-1);
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
	public static void deleteDatabases(File db){
		File d1=new File(db+".nhr");
		File d2=new File(db+".nin");
		File d3=new File(db+".nsq");
		if(d1.exists())d1.delete();
		if(d2.exists())d2.delete();
		if(d3.exists())d3.delete();

	}

}
