package blastTools;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;

public class PerformBlast {
	public static boolean blast(String blastpath,String builderpath,String program,double eValue,File output,File input,File database,boolean verbose,boolean createDB,boolean DNA,boolean Filter){
		try{
			File db=new File(database+".nhr");
			String dbType=DNA?"F":"T";
			if(!db.exists()||createDB){
				String blastcom=builderpath+" -p "+dbType+" -i "+database;
				Process p=Runtime.getRuntime().exec(blastcom);
				if(p.waitFor()!=0){
					System.err.println(blastcom);
					InputStream i=p.getErrorStream();
					int c=0;
					while((c=i.read())!=-1){
						System.err.print((char)c);
					}
					System.err.println("Formatdb was not successful!");
					return false;
				}
			}
			char filter=Filter?'T':'F';
			String blastcom=blastpath+" -p "+ program+" -m 8 -o "+output+" -i "+input+" -e "+eValue+" -d "+database+" -F "+filter+" -b 50000000";
			if(!verbose)System.out.println(blastcom);
			Process p=Runtime.getRuntime().exec(blastcom);
			if(p.waitFor()!=0){
				System.err.println(blastcom);
				InputStream i=p.getErrorStream();
				int c=0;
				while((c=i.read())!=-1){
					System.err.print((char)c);
				}
				System.err.println("Blast was not successful!");
				return false;
			}
		}catch(InterruptedException e){
			e.printStackTrace();
			return false;
		}catch(IOException e){
			e.printStackTrace();
			return false;
		}
		return true;
		
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
