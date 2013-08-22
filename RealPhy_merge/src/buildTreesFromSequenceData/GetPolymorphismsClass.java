package buildTreesFromSequenceData;

import java.io.*;
import java.util.*;
import java.util.Map.*;


import solexa.*;
import solexa.Arrays;
import util.*;


public abstract class GetPolymorphismsClass implements GetPolymorphisms,Serializable{
	/**
	 * 
	 */
	
	public final int BAM=2;
	public final int soap=0;
	public final int SAM=1;
	private static final long serialVersionUID = 1L;
	ArrayList<Genes> genePoly=new ArrayList<Genes>();
	HashMap<String /*strain*/,HashMap<String /*geneID*/,HashMap<Integer/*position*/,Polymorph/*poly*/>>> strains=new HashMap<String, HashMap<String,HashMap<Integer,Polymorph>>>();

	HashMap<String,HashMap<String,HashMap<Integer,Integer>[]>> baseNames=new HashMap<String, HashMap<String,HashMap<Integer,Integer>[]>>();
	File Reference=null;
	int flank=0;
	int quality=20;
	double polymorphismThreshold=0.95;
	double fractionCovThreshold=0.93;
	int perBaseCoverage=10;
	boolean genes=true;
	boolean printInvariant=false;
	//int covWindow=500;
	File outFolder;
	double gapThreshold=0;
	boolean allowGaps=true;
	String reference;
	private static final HashMap<Character,Boolean> ATGC=new HashMap<Character, Boolean>();

	static {
		ATGC.put('A',true);
		ATGC.put('T',true);
		ATGC.put('G',true);
		ATGC.put('C',true);
		ATGC.put('a',true);
		ATGC.put('t',true);
		ATGC.put('g',true);
		ATGC.put('c',true);
	}
	//int readLength;
	public static ArrayList<File> getSoapFiles(File inputFolder){
		ArrayList<File> soapFiles=new ArrayList<File>();
		File[] list=new File(inputFolder+"/soapout/").listFiles();
		for(int i=0;i<list.length;i++){
			if(list[i].isFile()&&list[i].getName().endsWith(".sop")){
				soapFiles.add(list[i]);
			}
		}
		return soapFiles;
	}
	/*
	private void print(TreeMap<Integer,Character> tm){
		Iterator <Entry<Integer,Character>> it=tm.entrySet().iterator();
		while(it.hasNext()){
			Entry<Integer,Character> e=it.next();
			System.out.println(e.getKey()+" "+e.getValue());
		}
	}*/
	
	  ArrayList<Integer> deleteUnsureRegions(PointSubstitutions pss,int geneNumber) {
		String id=genePoly.get(geneNumber).geneID;
		Iterator<Entry<Integer,Character>> it=genePoly.get(geneNumber).polymorphisms.entrySet().iterator();
		ArrayList<Integer> delete=new ArrayList<Integer>();
		double[] coverage=pss.getCoverage(id);
		Arrays geneBases=pss.getBases(id);
		int k=0;
		while(it.hasNext()){
			Entry<Integer,Character> e=it.next();
			int frame=genes?k%3:0;
			int pos=e.getKey();
			char base=e.getValue();
			
			double numMajorPoly=geneBases.numMajorPolymorphism(pos);
			double ratioPoly=numMajorPoly/(coverage[pos]*1.0);
			if((pos<coverage.length&&coverage[pos]<perBaseCoverage)||(pos>=coverage.length)||!ATGC.containsKey(base)||(ratioPoly < polymorphismThreshold && ratioPoly > (1- polymorphismThreshold ))){//delete whole triplets

				//if((pos<coverage.length&&coverage[pos]<perBaseCoverage)||(pos>=coverage.length)||!ATGC.containsKey(base)||!((numMajorPoly/(coverage[pos]*1.0) >= polymorphismThreshold ||(printInvariant&&numMajorPoly/(coverage[pos]*1.0) <= (1-polymorphismThreshold)))) ){	
				delete.add(pos-frame);
				if(genes){
					delete.add(pos-frame+1);
					delete.add(pos-frame+2);
				}
			}
			k++;
		}
		return delete;
	}
	
    private  int initGP(ArrayList<Fasta> fas,boolean addReference){
    	int sum=0;
		for(int i=0;i<fas.size();i++){
			String ident=fas.get(i).getIdent();
			String split[]=ident.split("\\s+");
			String id=split[0];
			String pos=split[1].split("\\.")[0];
			char orient=split[2].charAt(0);
			String contig=split[3];
            sum+=fas.get(i).getSequence().length();
            
			genePoly.add(new Genes(id,addReference?referenceTreeMap(fas.get(i).getSequence()):new TreeMap<Integer,Character>(),fas.get(i).getSequence().length()-2*flank,Integer.parseInt(pos),orient,contig));
		}
		
		return sum;
	}
    
    private TreeMap<Integer,Character> referenceTreeMap(String seq){
    	TreeMap<Integer,Character> tm=new TreeMap<Integer, Character>();
    	for(int i=0+flank;i<seq.length()-flank;i++){
    		tm.put(i+1-flank, seq.charAt(i));//changed
    	}
    	return tm;
    }
	public GetPolymorphismsClass(ArrayList<File> soapFiles,ArrayList<String> references,File ReferenceFas,int Flank,int Quality,double PolymorphismThreshold,double FractionCovThreshold,int PerBaseCoverage,boolean subInfo,boolean NoGenes,boolean PrintInvariant,File OutFolder,boolean addReference,double gapThreshold,String ref) {
		Reference=ReferenceFas;
		genes=!NoGenes;
		flank=Flank;
		quality=Quality;
		polymorphismThreshold=PolymorphismThreshold;
		fractionCovThreshold=FractionCovThreshold;
		perBaseCoverage=PerBaseCoverage;
		printInvariant=PrintInvariant;
		outFolder=OutFolder;
		this.gapThreshold=gapThreshold;
		calculatePolymorphisms(soapFiles,references,subInfo,addReference);
		reference=ref;

	}
    
	public  void writeGeneList(File out){
		try{

			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			bw.write("GeneID\tPos\tlength\tL/#P\tPoly\n");

			for(int i=0;i<genePoly.size();i++){
				int length=genePoly.get(i).length;
				double polyprop=genePoly.get(i).polymorphisms.size()/(1.0*length);
				bw.write(genePoly.get(i).geneID+"\t"+genePoly.get(i).pos+"\t"+length+"\t"+polyprop+"\t"+genePoly.get(i).polymorphisms.size()+"\n");
			}
			bw.close();
		}catch(IOException e){
			e.printStackTrace();
		}
	}

	public String getStrain(String name){
		String split[]=name.split("_");
		if(split[split.length-1].endsWith("fasta")){
			StringBuffer newName=new StringBuffer();
			newName.append(split[0]);
			for(int i=1;i<split.length-1;i++){
				newName.append("_"+split[i]);
			}
			return newName.toString();
		}
		return name;
	}
	
	private HashMap<String,Boolean> makeHash(ArrayList<String> list){
		HashMap<String,Boolean> hm=new HashMap<String, Boolean>();
		for(int i=0;i<list.size();i++){
			hm.put(list.get(i),true);
		}
		return hm;
	}
	
	void calculatePolymorphisms(ArrayList<File> alignmentFiles,ArrayList<String> references,boolean subInfo,boolean addReference){
		ArrayList<Fasta> ref=Fasta.readFasta(Reference);
        int refLength=initGP(ref,addReference);
        HashMap<String,Boolean> hmRef=makeHash(references);
        try{
        	BufferedWriter bw=new BufferedWriter(new FileWriter(outFolder+"/polies_and_orthologous_sites.txt"));
        	BufferedWriter bwCore=new BufferedWriter(new FileWriter(outFolder+"/coreGenomeSize.txt"));

        	System.out.println("Strain\t#mapped nucleotides/total possible|Number of polymorphic sites present in all strains analysed to that point\n");
        	bw.write("Name|mappedSites|polymorphisms|ratio\n");//|min % missed SNPs\n");
			TreeMap<Integer,File> sortedFiles=new TreeMap<Integer,File>();
			
        	for(int i=0;i<alignmentFiles.size();i++){
        		File alignmentFile=alignmentFiles.get(i);
        		int aligner=getAligner(alignmentFile);

        		String strainIntern=alignmentFile.getName().split("\\.")[0];
        		String strainExtern=getStrain(strainIntern);
        		boolean basenames=subInfo&&hmRef.containsKey(strainExtern);
        		PointSubstitutions pss=aligner==soap?new PointSubstitutionsSoap(Reference,flank,alignmentFile,quality,basenames):new PointSubstitutionsBAM(Reference,flank,alignmentFile,quality,basenames);
        		HashMap<String,HashMap<Integer,Polymorph>> strainPolies=adjustGP(pss,basenames,strainExtern);
        		strains.put(strainIntern,strainPolies);
        		int numberMappedSites=	pss.getMappedSites(perBaseCoverage);
        		int strainPolymorphisms=getStrainPolymorphisms(strainPolies);
        		double proportionMutations=strainPolymorphisms/(numberMappedSites*1.0);
        		bw.write(strainExtern+"|"+numberMappedSites+"|"+strainPolymorphisms+"|"+proportionMutations+"\n");//"|"+samplingBias+"\n");
        		

        		System.out.println(strainExtern+"\t"+numberMappedSites+"/"+refLength);//+(Runtime.getRuntime().totalMemory()/1048576)+"MB");
        		while(sortedFiles.containsKey(numberMappedSites)){
        			numberMappedSites+=1;
        		}
        		sortedFiles.put(numberMappedSites,alignmentFiles.get(i));
        		
        	}    
        	bw.close();
        	System.out.println("Overall quality check. Strains sorted by divergence from reference strain...");
        	Iterator<Entry<Integer,File>> it=sortedFiles.descendingMap().entrySet().iterator();
        	while(it.hasNext()){
        		Entry<Integer,File> e=it.next();
        		File alignmentFile=e.getValue();
        		String strainIntern=alignmentFile.getName().split("\\.")[0];
        		String strainExtern=getStrain(strainIntern);
        		int aligner=getAligner(alignmentFile);
        		boolean basenames=subInfo&&hmRef.containsKey(strainExtern);
        		PointSubstitutions pss=aligner==soap?new PointSubstitutionsSoap(Reference,flank,alignmentFile,quality,basenames):new PointSubstitutionsBAM(Reference,flank,alignmentFile,quality,basenames);

        		checkCoverageGP(pss,strainIntern);
        		if(basenames){
        			baseNames.put(strainIntern, pss.getBaseNames());
        		}
        		int numberPolies=getNumberPolies();
        		int mappedSites=pss.getMappedSites(perBaseCoverage);
        		System.out.println(strainExtern+"\t"+mappedSites+"/"+refLength+"|"+numberPolies);//"|"+(Runtime.getRuntime().totalMemory()/1048576)+"MB");
        		bwCore.write(strainExtern+"\t"+mappedSites+"/"+refLength+"|"+numberPolies+"\n");
        	}
        	writeCoverage(new File(outFolder+"/coverage.txt"));
        	bwCore.close();
        }catch(IOException e){
        	e.printStackTrace();
        }

	}
	 
	private int getAligner(File alignmentFile){
		if(alignmentFile.toString().endsWith(".sop")){
			return soap;
		}else if (alignmentFile.toString().endsWith(".sam")){
			return SAM;
		}else{
			return BAM;
		}
	}
	
	 private static int getStrainPolymorphisms(HashMap<String,HashMap<Integer,Polymorph>> hm){
		 int count=0;
		 Iterator<Entry<String,HashMap<Integer,Polymorph>>> it=hm.entrySet().iterator();
		 while(it.hasNext()){
			 Entry<String,HashMap<Integer,Polymorph>> e=it.next();
			 count+=e.getValue().size();
		 }
		 return count;
	 }

	 /*needs changing!!!
	  * 
	  */
	 
	 public int getNumberPolies(){
		 int sum=0;
		 int size=genePoly.size();
		 for(int i=0;i<size;i++){
			 sum+=genePoly.get(i).polymorphisms.size();
		 }
		 return sum;
	 }
	 public void writeCoverage(File out){
		 try{
			 BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			 int size=genePoly.size();
			 for(int i=0;i<size;i++){
				 Genes gene=genePoly.get(i);
				 int geneLength=gene.length;
				 for(int j=1;j<geneLength+1;j++){
					 if(gene.polymorphisms.containsKey(j)){
						 
						 bw.write(gene.geneID+"\t1\n");
					 }else{
						 bw.write(gene.geneID+"\t0\n");
					 }
				 }
			 }
			 bw.close();
		 }catch(IOException e){
			 e.printStackTrace();
		 }
	 }
	 
	public String[] getStrain(String[] list){
		String[] newList=new String[list.length];
		for(int i=0;i<newList.length;i++){
			newList[i]=getStrain(list[i]);
		}
		return newList;
	}
	
	private int getReferencePos(ArrayList<String> sorted){
		for(int i=0;i<sorted.size();i++){
			//System.out.println(sorted.get(i));
			if(sorted.get(i).equals(reference)){
				return i;
			}

		}
		System.err.println("Could not find reference "+reference);
		throw new RuntimeException();
	}
	
	private int getReferencePos(ArrayList<String> sorted,String[] sortedInt){
		int refNum=0;
		for(int i=0;i<sorted.size();i++){
			//System.out.println(sorted.get(i));
			if(sorted.get(i).equals(reference)){
				return refNum;
			}
			if(baseNames.containsKey(sortedInt[i])){
				refNum++;
			}
		}
		System.err.println("Could not find reference "+reference);
		throw new RuntimeException();
	}
	
	public Clashes calculateColumns(){
		
		
		String[] strainListIntern=strains.keySet().toArray(new String[0]);
		TreeMap<String,Boolean> sortedStrains=fillTreeMap(strainListIntern);
		String[] sorted=sortedStrains.keySet().toArray(new String[0]);
		ArrayList<String> sortedEx=toExtern(sorted);
		int ref=getReferencePos(sortedEx,sorted);
	//	int refSorted=getReferencePos(sortedEx);
		Clashes columns=new Clashes(sortedEx);
		for(int i=0;i<genePoly.size();i++){
			Genes gene=genePoly.get(i);
			Iterator<Entry<Integer,Character>> it2= gene.polymorphisms.entrySet().iterator();
			while(it2.hasNext()){
				
				Entry<Integer,Character> e2=it2.next();
				int pos=e2.getKey();
				char origBase=e2.getValue();
				StringBuffer baseColumn=new StringBuffer();
				ArrayList<Integer> queryIDColumn=new ArrayList<Integer>();
				boolean addition=false;
				for(int j=0;j<sorted.length;j++){
					String strainEx=getStrain(sorted[j]);
					HashMap<Integer,Polymorph> polies=strains.get(sorted[j]).get(gene.geneID);
					char base=polies.containsKey(pos)?polies.get(pos).base:origBase;
					baseColumn.append(base);
					
					//if(polies.get(pos).name.equals("S21_39129")
					
					
					//id it is not one of the reference strains then skip this step
					if(!baseNames.containsKey(sorted[j]))continue;
					HashMap<Integer,Integer>[] GeneBaseNames=baseNames.get(sorted[j]).get(gene.geneID);

					int queryID;
					if(polies.containsKey(pos)){
						if(polies.get(pos).name==null){
							System.err.println("something went wrong! polies get pos");
							System.err.println("pos: "+pos+" strain: "+strainEx);
							System.exit(-1);
						}
						queryID=polies.get(pos).name;
						addition=true;
					}else{
						if(GeneBaseNames[pos]==null){
							System.err.println("something went wrong!");
							System.exit(-1);
						}
						queryID=getMajority(GeneBaseNames[pos]);
					}
					
					queryIDColumn.add(queryID);
				}
				if(printInvariant||addition){
					columns.addColumn(queryIDColumn, baseColumn,(byte)ref);
				}
			}
		}
		return columns;

	}
	int getMajority(HashMap<Integer, Integer> hashMap){
		Iterator<Entry<Integer,Integer>> it=hashMap.entrySet().iterator();
		int max=0;
		int maxItem=-1;
		while(it.hasNext()){
			Entry<Integer,Integer> e=it.next();
			if(e.getValue()>max){
				max=e.getValue();
				maxItem=e.getKey();
			}
			
		}
		return maxItem;
		
	}
	
	private ArrayList<String> toExtern(String[] internIDs){
		ArrayList<String> strains=new ArrayList<String>();
		for(int i=0;i<internIDs.length;i++){
			strains.add(getStrain(internIDs[i]));
		}
		return strains;
	}
	
/*	private void print(HashMap<String,HashMap<String,String[]>> hm){
		Iterator<Entry<String,HashMap<String,String[]>>> it=hm.entrySet().iterator();
		while(it.hasNext()){
			System.out.println(it.next().getKey());
		}
	}*/

	private TreeMap<String,Boolean> fillTreeMap(String[] strains){
		TreeMap<String,Boolean> sorted=new TreeMap<String, Boolean>();
		for(int i=0;i<strains.length;i++){
			sorted.put(strains[i],true);
		}
		return sorted;
	}

}
 