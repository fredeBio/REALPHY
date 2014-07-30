package buildTreesFromSequenceData;

import java.io.*;
import java.util.*;

import util.*;
import blastTools.*;

public class DetermineCoreGenome {
	public static void main(String args[]){
		runDetermineCoreGenome(args);
	}

	private static File getGBKFile(String prefix,String[] ext){
		for(int i=0;i<ext.length;i++){
			File gbk=new File(prefix+"."+ext[i]);
			if(gbk.exists()){
				return gbk;
			}
		}
		return null;
	}
	
	public static ArrayList<File> runDetermineCoreGenome(String args[]){
		File inFolder=new File(args[0]);
		File outFolder=new File(args[1]);
		double perBaseScore=Double.parseDouble(args[2]);
		int flank=Integer.parseInt(args[3]);
		String program =args[4];
		int numberGenomes=Integer.parseInt(args[5]);
		int numberRepeats=Integer.parseInt(args[6]);
		String BLASTPATH=args[7];
		String BLASTBUILDER=args[8];
		String ref=args[9];
		ArrayList<File> cores=new ArrayList<File>();
		for(int l=0;l<numberRepeats;l++){
			try{
				//HashMap<String,Boolean> genesHM=new HashMap<String, Boolean>();
				int k=0;


				boolean DNA=program.equalsIgnoreCase("blastn");
				File[] list=inFolder.listFiles();			
				ArrayList<String> ids=chooseGenomes(numberGenomes,list,ref);



				String coreId=nameCore(ids,outFolder);
				File core=new File(coreId+".fas");
				File coreFlank=new File(coreId+"Flank"+flank+".fas");
				cores.add(coreFlank);
				HashMap<String,Fasta> coreHash=new HashMap<String, Fasta>();
				HashMap<String,Fasta> coreHashFlank=new HashMap<String, Fasta>();
				for(int i=0;i<ids.size();i++){
					String id=ids.get(i);
					System.out.println(id);

					File gbk=getGBKFile(inFolder+"/"+id,RealPhy.gbkExt);
					if(gbk==null){
						System.err.println("Cannot find genbank file. Probably an error in the code!");
						System.exit(-1);
					}
					//HashMap<String,String> fas=Fasta.fasToHash(Fasta.readFasta(new File(inFolder+"/"+id+".fas")));
					k++;
					if(k==1){
							coreHashFlank=ConvertGBKToFasta.writeCDS(gbk,  coreFlank, flank, DNA);
							coreHash=ConvertGBKToFasta.writeCDS(gbk,  core, 0, DNA);
					}else{
						File gFile=new File(outFolder+"/gFile.fas");
						gFile.deleteOnExit();
						ConvertGBKToFasta.writeCDS(gbk,gFile, 0, DNA);
						File blastout=new File(outFolder+"/blast.out");
						blastout.deleteOnExit();
						PerformBlast.blast(BLASTPATH,BLASTBUILDER,program, 1e-10, blastout, core, gFile, true, true,DNA,false);
						File tempNoFlank=new File(outFolder+"/tempNoFlank.fas");
						File tempFlank=new File(outFolder+"/temp.fas");
						BufferedWriter bwFlank=new BufferedWriter(new FileWriter(tempFlank));
						BufferedWriter bwNoFlank=new BufferedWriter(new FileWriter(tempNoFlank));
						ReadBlast rb=new ReadBlast(blastout);
						HashMap<String,Fasta> tempHashFlank=new HashMap<String, Fasta>();
						HashMap<String,Fasta> tempHash=new HashMap<String, Fasta>();
						String old="";
						for(int j=0;j<rb.getQuery().size();j++){
							if(!old.equals(rb.getQuery().get(j))){
								old=rb.getQuery().get(j);
								String gene=coreHash.get(old).getSequence();
								if(rb.getScore().get(j)/gene.length()>perBaseScore){
									bwFlank.write(">"+coreHash.get(old).getIdent()+"\n"+coreHashFlank.get(old).getSequence()+"\n");
									bwNoFlank.write(">"+coreHash.get(old).getIdent()+"\n"+coreHash.get(old).getSequence()+"\n");
									tempHash.put(old, coreHash.get(old));
									tempHashFlank.put(old, coreHashFlank.get(old));
								}
							}
						}
						bwFlank.close();
						bwNoFlank.close();
						coreHash=tempHash;
						coreHashFlank=tempHashFlank;
						tempFlank.renameTo(coreFlank);
						tempNoFlank.renameTo(core);
					}

				}
			}catch(IOException e){
				e.printStackTrace();
				System.exit(-1);
			}
		}
		return cores;
	}
	
	public static String nameCore(ArrayList<String> ids,File outFolder){
		StringBuffer sb=new StringBuffer(outFolder+"/core+"+ids.get(0));
//		for(int i=1;i<ids.size();i++){
//			sb.append("+"+ids.get(i));
//		}
		return sb.toString();
	}
	
	public static ArrayList<String> chooseGenomes(int genomes,File[] list,String ref){
		
		ArrayList<String> ids=new ArrayList<String>();
		ArrayList<File> alist=new ArrayList<File>();
		int refNumber=-1;
		for(int i=0;i<list.length;i++){
			String name=list[i].getName();
			if(RealPhy.hasExtension(name,RealPhy.gbkExt)||RealPhy.hasExtension(name,RealPhy.fasExt)){
				alist.add(list[i]);
				if(ref.equals(list[i].getName().split("\\.")[0])){
					refNumber=alist.size()-1;
				}
			}
		}
		int start=0;
		if(ref.length()>0){
			ids.add(ref);
			alist.remove(refNumber);
			start=1;
		}
		for(int i=start;alist.size()>0&&i<genomes;i++){
			int rand=(int)(Math.random()*alist.size());
			File current=alist.get(rand);
			String id=current.getName().split("\\.")[0];
			ids.add(id);
			alist.remove(rand);
		}
		return ids;
	}
}
