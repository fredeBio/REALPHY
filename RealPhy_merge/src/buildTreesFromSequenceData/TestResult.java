package buildTreesFromSequenceData;

import java.io.*;
import java.util.*;

import util.*;
import util.phylogenetics.*;

public class TestResult {
	public static void main(String args[]){
		File alignment=new File(args[0]);
		File detailsFolder=new File(args[1]);
		String ref=args[2];
		File out=new File(args[3]);
		boolean genes=Boolean.parseBoolean(args[4]);
		if(genes){
			ArrayList<ArrayList<Fasta>> seqs=getSequences(alignment.getParentFile());
			alignment=makeTrueAlignment(alignment,seqs);
		}
		TestResult tr=new TestResult(detailsFolder, alignment, ref,out);
	}
	
	/*private String[] combine(String[] list1,String[] list2){
		String combined[]=new String[list1.length+list2.length];
		for(int i=0;i<list1.length;i++){
			combined[i]=list1[i];
		}
		for(int i=list1.length;i<combined.length;i++){
			combined[i]=list2[i-list1.length];
		}
		return combined;
	}*/
	
	private static  ArrayList<ArrayList<Fasta>> getSequences(File seqfolder){
		File[] list=seqfolder.listFiles();
		ArrayList<ArrayList<Fasta>> seqs=new ArrayList<ArrayList<Fasta>>();
		for(int i=0;i<list.length;i++){
			String name=list[i].getName();
				if(RealPhy.hasExtension(name, RealPhy.gbkExt)){
					ReadGenbank rgb=new ReadGenbank(list[i]);
					ArrayList<Fasta> fas=rgb.getSequence();
					seqs.add(fas);
				}else if(RealPhy.hasExtension(name, RealPhy.fasExt)){
					ArrayList<Fasta> fas=Fasta.readFasta(list[i]);
					seqs.add(fas);
				}
		}
		return seqs;
		
	}
	
	/**
	 * TODO! ONLY WORKS FOR GENEDEFS IN A SINGLE FASTA SEQUENCE (IE NO CONTIGS ETC)
	 * @param gbk
	 * @param seqs
	 * @return
	 */
	
	private static File makeTrueAlignment(File gbk,ArrayList<ArrayList<Fasta>> seqs){
		ReadGenbank rgb=new ReadGenbank(gbk,"CDS");
		ArrayList<Info> regions=rgb.getFirstSequenceFeatures();
		ArrayList<Fasta> trueAl=new ArrayList<Fasta>();
		File out=new File(gbk.getParent()+"/trueAlignment.fas_temp");
		for(int i=0;i<seqs.size();i++){
			String id=seqs.get(i).get(0).getIdent();
			StringBuffer temp=new StringBuffer();
			String seq=seqs.get(i).get(0).getSequence();
			for(int j=0;j<regions.size();j++){
				int start=regions.get(j).getStart();
				int end=regions.get(j).getEnd();
				temp.append(seq.substring(start,end));
			}
			trueAl.add(new Fasta(id,temp.toString()));
		}
		Fasta.write(trueAl, out);
		return out;
	}
	
	
	Alignment algTrue;
	Alignment algRP;
	ArrayList<Integer> cov=new ArrayList<Integer>();
	ArrayList<Integer> alignmentPos=new ArrayList<Integer>();
	HashMap<String,Integer> transRP;
	HashMap<String,HashMap<Integer,Character>> details=new HashMap<String, HashMap<Integer,Character>>();
	int referencePos;
	File outFolder;
	public TestResult(File detailsFolder,File alignment,String ref,File outFolder){
			algTrue=new Alignment(Fasta.readFasta(alignment));
			algRP=new Alignment(Fasta.readFasta(new File(detailsFolder+"/polymorphisms_move.fas")));
			referencePos=getReference(ref);
			transRP=trans(algRP.getIdents());
			this.outFolder=outFolder;
			File coverageFile=new File(detailsFolder+"/coverage.txt");
			readCoverage(coverageFile);
			readDetails(detailsFolder);
			testRealPhy();
			makeAlignment();
			compareAlignments();
		
	}
	
	private void makeAlignment(){
			File out=new File(outFolder+"/newAlignment.fas");
			Alignment alg=new Alignment();
			int j=0;
			for(int i=0;i<cov.size();i++){
				if(cov.get(i)==1){
					alg.addColumn(algRP.getColumn(j));
					j++;
				}else{
					alg.addColumn("---");
				}
			}
			Fasta.write(alg.toFasta(), out);
	}
	
	public void printIdents(Alignment alg){
		ArrayList<String> idents=alg.getIdents();
		for(int i=0;i<idents.size();i++){
			System.out.println(idents.get(i));
		}
	}
	public void compareAlignments(){
		int sum=0;

		for(int i=0;i<alignmentPos.size();i++){
			//System.out.println(algRP.getColumn(i));
			//System.out.println(algTrue.getColumn(alignmentPos.get(i)+1));
			if(!equals(i,alignmentPos.get(i))){
				sum++;
				System.out.println("Alignment column "+i+" in RP alignment does not equal true alignment at pos "+alignmentPos.get(i)+".");
				//System.exit(-1);
			}
		}
		System.out.println(sum);
	}
	public void readDetails(File detailsFolder){
		ArrayList<String> idents=algTrue.getIdents();
		for(int i=0;i<idents.size();i++){
			File d=new File(detailsFolder+"/"+idents.get(i)+"details.txt");
			details.put(idents.get(i),readDetail(d));
			
		}
 	}
	
	public HashMap<Integer,Character> readDetail(File detail){
		HashMap<Integer,Character> hm=new HashMap<Integer, Character>();

		try{
			BufferedReader br=new BufferedReader(new FileReader(detail));
			String line;
			int i=0;
			while((line=br.readLine())!=null){
				i++;
				if(i==1)continue;
				String split[]=line.split("\\s+");
				int truePos=Integer.parseInt(split[5])-1;
				hm.put(truePos, split[7].toUpperCase().charAt(0));
			}
			br.close();
		}catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
		return hm;
	}
	
	private int getReference(String ref){
		ArrayList<String> idents=algTrue.getIdents();
		for(int i=0;i<idents.size();i++){
			if(ref.equals(idents.get(i))){
				return i;
			}
		}
		return -1;
	}
	public void testRealPhy(){
		if(referencePos==-1){
			System.err.println("Reference not found in alignment.");
			System.exit(-1);
		}
		int sum=0;
		for(int i=0;i<cov.size();i++){
			if(cov.get(i)==1){
				sum++;
				boolean invarAlg=isInvariantAlg(i);
				boolean invarDet=isInvariantDetails(i);
				if(invarAlg!=invarDet){
					System.err.println("Problem with alignment position "+i+". Alignment site is "+(invarAlg?"invariant":"variable")+" and site in details file is "+(invarDet?"invariant.":"variable."));
				}else{
					if(!invarAlg){
						if(!isEqualVariable(i)){
							System.err.println("Called polymorphisms differ between correct alignment and details file at position: "+i);
						}
					}
				}
			}
		}
		System.out.println(sum);
	}
	
	private boolean isEqualVariable(int pos){
		ArrayList<String> idents=algTrue.getIdents();
		String site=algTrue.getColumn(pos).toUpperCase();
		char refC=site.charAt(referencePos);
		for(int i=0;i<idents.size();i++){
			if(i==referencePos)continue;
			char cAlg=site.charAt(i);
			char cDet;
			if(details.get(idents.get(i)).containsKey(pos)){
				cDet=details.get(idents.get(i)).get(pos);
			}else{
				cDet=refC;
			}
			if(cDet!=cAlg){
				return false;
			}
		}
		return true;
	}
	private HashMap<String, Integer> trans(ArrayList<String> id){
		HashMap<String, Integer> hm=new HashMap<String, Integer>();
		for(int i=0;i<id.size();i++){
			hm.put(id.get(i), i);
		}
		return hm;
	}
	private boolean equals(int posRP,int posTrue){
		ArrayList<String> idents=algTrue.getIdents();
		String siteTrue=algTrue.getColumn(posTrue).toUpperCase();
		String siteRP=algRP.getColumn(posRP).toUpperCase();
		//System.out.println(siteTrue+" "+siteRP);
		for(int i=0;i<idents.size();i++){
			char cTrue=siteTrue.charAt(i);
			char cRP=siteRP.charAt(transRP.get(idents.get(i)));
			
			if(cRP!=cTrue){
				return false;
			}
		}
		return true;
	}
	
	private boolean isInvariantDetails(int pos){
		ArrayList<String> idents=algTrue.getIdents();
		for(int i=0;i<idents.size();i++){
			if(details.get(idents.get(i)).containsKey(pos)){
				return false;
			}
		}
		return true;
	}
	
	private boolean isInvariantAlg(int pos){
		String site=algTrue.getColumn(pos);
		char c=site.charAt(0);
		for(int i=1;i<site.length();i++){
			if(c!=site.charAt(i)){
				return false;
			}
		}
		return true;
	}
	
	public void readCoverage(	File coverage){
		try{
			
			BufferedReader br=new BufferedReader(new FileReader(coverage));
			String line;
			int i=0;
			while((line=br.readLine())!=null){
				String split[]=line.split("\\s+");
				int covered=Integer.parseInt(split[1]);
				cov.add(covered);
				if(covered==1){
					alignmentPos.add(i);
				}
				i++;
			}
			
			br.close();
		}catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
}
