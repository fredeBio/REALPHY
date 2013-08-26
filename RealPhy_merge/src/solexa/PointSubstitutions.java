package solexa;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map.Entry;

import util.DNAmanipulations;
import util.Fasta;
import util.Info;

public abstract class PointSubstitutions {
	
	public PointSubstitutions(File RefSeq,int flank,File AlignmentFile,int quality,boolean subInfo){
		this(RefSeq,flank,AlignmentFile,quality,1,subInfo);
	}
	public PointSubstitutions(File RefSeq,int flank,File AlignmentFile,int quality,int fold,boolean subInfo){
		//readSeq(RefSeq, flank);
		fasta=readFasta(RefSeq,flank);
		coverage = initCoverage();
		coveragePos = initCoverage();
		coverageNeg = initCoverage();
		baseNames=initNames();
		substitutions = initCoverage();
		bases = initBases();
		alignmentFile=AlignmentFile;
		this.subInfo=subInfo;
		read( quality,flank,fold);
		for(int i=0;i<alName.size();i++){
			bases.get(alName.get(i)).setCoverage(coverage.get(alName.get(i)));
		}
		setAvgCoverage();
		setfractionCoverage();
	}
	//needs to call setCoverage and setSubstitutions
	abstract void read(int quality,int flank,int fold);
	
	
	
	ArrayList<Integer> alLength = new ArrayList<Integer>();
	ArrayList<String> alName = new ArrayList<String>();
	HashMap<String,HashMap<String,Integer>> substitutionInfo=new HashMap<String, HashMap<String,Integer>>();
	HashMap<String,Arrays> bases;
	HashMap<String,double[]> substitutions;
	HashMap<String,double[]> coverage;
	HashMap<String,double[]> coveragePos;
	HashMap<String,double[]> coverageNeg;
	HashMap<String,HashMap<Integer,Integer>[]> baseNames;

	HashMap<String,String> fasta;
	HashMap<String,Double> avgCov=new HashMap<String, Double>();
	HashMap<String,Double> fractionCov=new HashMap<String,Double>();
	File alignmentFile;
	boolean subInfo=false;
	
	public static ArrayList<Integer> getMismatchesRef(String mismatch){
		String pos[]=mismatch.split("[^0-9]+");
		String edit[]=mismatch.split("[0-9]+");
		ArrayList<Integer> substitutionPos=new ArrayList<Integer>();
		int start=0;
		for(int i=1;i<pos.length;i++){
			start=start+Integer.parseInt(pos[i-1])+1;
			if(edit[i].charAt(0)=='^'){
			}else{
				substitutionPos.add(start);
				
			}
		}
		return substitutionPos;
	}
	
	public static ArrayList<Integer> getMismatchesRead(String mismatch){
		String pos[]=mismatch.split("[^0-9]+");
		String edit[]=mismatch.split("[0-9]+");
		ArrayList<Integer> substitutionPos=new ArrayList<Integer>();
		int start=0;
		for(int i=1;i<pos.length;i++){
			start=start+Integer.parseInt(pos[i-1])+1;
			if(edit[i].charAt(0)=='^'){
				start--;
			}else{

				substitutionPos.add(start);
				
			}
		}
		return substitutionPos;
	}
	
	
	public static ArrayList<Integer> getMismatchesRead(String mismatch,ArrayList<Integer> refDeletion){
		String pos[]=mismatch.split("[^0-9]+");
		String edit[]=mismatch.split("[0-9]+");
		ArrayList<Integer> substitutionPos=new ArrayList<Integer>();
		int start=0;
		for(int i=1;i<pos.length;i++){
			start=start+Integer.parseInt(pos[i-1])+1;
			if(edit[i].charAt(0)=='^'){
				start--;
			}else{
				//System.out.println(start);

				substitutionPos.add(start);
				
			}
		}
		int refDelSize=refDeletion.size();
		int subPosSize=substitutionPos.size();
		int j=0;
		for(int i=0;i<refDelSize&&j<subPosSize;i++){
			if(refDeletion.get(i)<substitutionPos.get(j)){
				for(int k=j;k<subPosSize;k++){
					substitutionPos.set(k,substitutionPos.get(k)+1);
				}
			}else{
				j++;
				i--;
			}
		}
		return substitutionPos;
	}

	public static ArrayList<Integer> getRefDeletions(String cigar){ //reads cigar string and returns matching regions
		String[] distances=cigar.split("[^0-9]+");
		String[] editCodes=cigar.split("[0-9]+");
		ArrayList<Integer> positions=new ArrayList<Integer>();
		int start=0;
		int length=distances.length;
		for(int i=1;i<length+1;i++){
			//System.out.println(distances[i-1]+" "+editCodes[i]);
			if(editCodes[i].equals("D")||editCodes[i].equals("P")){
				start=start+Integer.parseInt(distances[i-1]);
			}else if(editCodes[i].equals("I")){
				int deletions=Integer.parseInt(distances[i-1]);
				for(int j=0;j<deletions;j++){
					start++;
					positions.add(start);
				}
			}else {
				start=start+Integer.parseInt(distances[i-1]);
			}
		}
		return positions;
	}

	
	public static void getMatchRegionsCigar(String cigar,StringBuffer read,StringBuffer ref){ //reads cigar string and returns matching regions
		String[] distances=cigar.split("[^0-9]+");
		String[] editCodes=cigar.split("[0-9]+");

		
		for(int i=1;i<distances.length+1;i++){
			int dist=Integer.parseInt(distances[i-1]);
			//match/substitution: move along in read and reference
			if(editCodes[i].equals("M")){
				for(int j=0;j<dist;j++){
					read.append('m');
					ref.append('m');
				}
				//if there is a deletion in the read, then only move along in the reference!
			}else if(editCodes[i].equals("D")){
				for(int j=0;j<dist;j++){
					read.append('g');
					ref.append('m');
				}
			//insertion in read or deletion in read overlapping another read insertion then move along in read coverage	
			}else if(editCodes[i].equals("I")){
				for(int j=0;j<dist;j++){
					read.append('m');
					ref.append('g');
				}
			}else{
			// if there is a deletion both in the read and the reference do nothing
			}
		}

	}
	
	
	public HashMap<String,Arrays> getBases(){
		return bases;
	}
	public Arrays getBases(String id){
		return bases.get(id);
	}
	public HashMap<String,double[]> getCoverage(){
		return coverage;
	}
	public HashMap<String,double[]> getCoveragePos(){
		return coveragePos;
	}
	public HashMap<String,double[]> getCoverageNeg(){
		return coverageNeg;
	}
	public double[] getCoverage(String id){
		return coverage.get(id);
	}
	public double[] getCoveragePos(String id){
		return coveragePos.get(id);
	}
	public double[] getCoverageNeg(String id){
		return coverageNeg.get(id);
	}
	public String getGene(String id){
		return fasta.get(id);
	}
	public HashMap<String,HashMap<Integer,Integer>[]> getBaseNames(){
		return baseNames;
	}
	public HashMap<Integer,Integer> getBaseName(String id,int pos){
		return baseNames.get(id)[pos];
	}
	
	public boolean checkfractionCoverage(String id,double t){
		
		return fractionCov.get(id)>t;
	}
	public double getfractionCoverage(String id){
		
		return fractionCov.get(id);
	}
	public boolean checkAvgCoverage(String id,double t){
		return avgCov.get(id)>t;
	}
	
	public void setAvgCoverage(){
		for(int i=0;i<alName.size();i++){
			int sum=0;
			for(int j=0;j<alLength.get(i);j++){
				sum+=coverage.get(alName.get(i))[j];
			}
			avgCov.put(alName.get(i), (sum/(alLength.get(i)*1.0)));
		}
		
	}
	
	public void setfractionCoverage(){
		for(int i=0;i<alName.size();i++){
			int sum=0;
			for(int j=0;j<alLength.get(i);j++){
				if(coverage.get(alName.get(i))[j]>0){
					sum++;
				}
			}
			fractionCov.put(alName.get(i), (sum/(alLength.get(i)*1.0)));
		}
		
	}
	
	public HashMap<String,String> readFasta(File in,int flank){
 		ArrayList<Fasta> fas=Fasta.readFasta(in);
		HashMap<String,String> hm=Fasta.fasToHash(fas,false);
		int size=fas.size();
		for(int i=0;i<size;i++){
			Fasta fasta=fas.get(i);
			alName.add(fasta.getIdent().split("\\s+")[0]);
			alLength.add(fasta.getSequence().length()-2*flank);
		}
		return hm;
	}
	
	public  HashMap<String,Arrays> initBases(){
		HashMap<String,Arrays> bases=new HashMap<String, Arrays>();
		for(int i=0;i<alLength.size();i++){
			bases.put(alName.get(i),new Arrays(alLength.get(i)+1));
		}
		return bases;
	}

	public  HashMap<String, HashMap<Integer, Integer>[]> initNames(){
		HashMap<String,HashMap<Integer,Integer>[]> baseNames=new HashMap<String, HashMap<Integer,Integer>[]>();
		for(int i=0;i<alName.size();i++){
			baseNames.put(alName.get(i), new HashMap[alLength.get(i)+1]);
		}
		return baseNames;
	}
	public HashMap<String, double[]> initCoverage(){
		HashMap<String,double[]> coverage=new HashMap<String, double[]>();
		for(int i=0;i<alName.size();i++){
			coverage.put(alName.get(i), new double[alLength.get(i)+1]);
		}
		return coverage;
	}
	
	private double initWindowCount(ArrayList<Double> countArray,int stepsize,int covWindowSize){
		double count=0;
		for(int i=0;i<covWindowSize/stepsize&&i<countArray.size();i++){
			count+=countArray.get(i);
		}
		return count;
	}
	
	private void addArea(int start,int end,ArrayList<Info> areas){
		if(areas.size()==0){
			Info temp=new Info(start,end,"");
			areas.add(temp);
		}else{
			Info previous=areas.get(areas.size()-1);
			if(previous.getEnd()>=start){
				Info temp=new Info(previous.getStart(),end,"");
				areas.remove(areas.size()-1);
				areas.add(temp);
			}else{
				Info temp=new Info(start,end,"");
				areas.add(temp);
			}
		}
	}
	
	private ArrayList<Double> getCountArray(int stepSize,double[] cov){
		ArrayList<Double> countArray=new ArrayList<Double>();
		double count=0;
		for(int i=0;i<cov.length;i++){
			count+=cov[i]>0?1:0;
			if((i+1)%10==0){
				countArray.add(count);
				count=0;
			}
		}
		return countArray;
	}
	
	public void writeArray( double threshold, File refseq, File outDir,int flank) {
		try {
			
			BufferedWriter bwCov=new BufferedWriter(new FileWriter(new File(outDir+"/coverage.out")));
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File(
					outDir + "/" + refseq.getName() + ".out")));
			int occ=0;
			int changes=0;
			bw.write("Pos\tcov\torigF\tF\tA\tT\tG\tC\tR\tA\tT\tG\tC\tGene\n");
			for(int i=0;i<alName.size();i++){
				String gene=fasta.get(alName.get(i));
				Arrays geneBases=bases.get(alName.get(i));
				double[] geneCoverage=coverage.get(alName.get(i));
				for(int j=0;j<geneCoverage.length;j++){
					System.out.println("cov "+geneCoverage[j]+"|");
					System.out.println("pol "+geneBases.numBases(j)+"|");

				}
				System.out.println();

				int cov=0;
				//System.out.println(geneBases.length);
				System.out.println(geneBases.length);
				System.out.println(threshold);

				for (int j = 1; j < geneBases.length-1-2*flank; j++) {
					cov+=geneCoverage[j];
					//System.out.println(j);
					if (geneBases.numBases(j)/(geneCoverage[j]*1.0) > threshold && geneCoverage[j]>0){//  && cov[i]>bases.numBases(i)+0.1*cov[i]) {//exclude position with less than threshold changes and actual substitutions
					//if  (cov[i]<=bases.numBases(i)+1) {	//only substitutions
					occ++;
					//minus 1
					/*bw.write( alName.get(i) + "\t" + (j-1) + "\t" + geneCoverage[j-1] + "\t" + gene.charAt(j-2)
								+ "\t\t" + geneBases.get(j-1, "AF") + "\t"
								+ geneBases.get(j-1, "TF") + "\t" + geneBases.get(j-1, "GF")
								+ "\t" + geneBases.get(j-1, "CF") + "\t\t"
								+ geneBases.get(j-1, "AR") + "\t" + geneBases.get(j-1, "TR")
								+ "\t" + geneBases.get(j-1, "GR") + "\t"
								+ geneBases.get(j-1, "CR") +"\n");*/
					//plus 1
					/*bw.write(alName.get(i)+"\t"+j+1 + "\t" + geneCoverage[j+1] + "\t" + gene.charAt(j)
								+ "\t\t" + geneBases.get(j+1, "AF") + "\t"
								+ geneBases.get(j+1, "TF") + "\t" + geneBases.get(j+1, "GF")
								+ "\t" + geneBases.get(j+1, "CF") + "\t\t"
								+ geneBases.get(j+1, "AR") + "\t" + geneBases.get(j+1, "TR")
								+ "\t" + geneBases.get(j+1, "GR") + "\t"
								+ geneBases.get(j+1, "CR") +"\n");*/
					//actual base
					bw.write(geneBases.numBases(j)+"\t"+alName.get(i)+"\t"+j + "\t" + geneCoverage[j] + "\t" + gene.charAt(j-1+flank)
							+ "\t\t" + geneBases.get(j, "AF") + "\t"
							+ geneBases.get(j, "TF") + "\t" + geneBases.get(j, "GF")
							+ "\t" + geneBases.get(j, "CF") + "\t\t"
							+ geneBases.get(j, "AR") + "\t" + geneBases.get(j, "TR")
							+ "\t" + geneBases.get(j, "GR") + "\t"
							+ geneBases.get(j, "CR") +"\t");

//					int frame=(j-1)%3;
//					String triplet=gene.substring(j-frame-1+flank,j+3-frame-1+flank);
//					ArrayList<String> alternate=getAlt(triplet,j,Math.abs(frame),geneBases,false,3);
//					String origAA=DNAmanipulations.translate(triplet,DNAmanipulations.code());
//					bw.write("\toriginal codon: "+triplet+"="+origAA+" alternatives: ");
//
//					boolean altern=false;
//					for(int k=0;k<alternate.size();k++){
//						String alt=DNAmanipulations.translate(alternate.get(k),DNAmanipulations.code());
//						if (!alt.equals(origAA)){
//							altern=true;
//							bw.write(alternate.get(k)+"="+alt+" ");
//						}
//					}
//					if(altern)changes++;



					bw.write("\n");
				}
				}
				bwCov.write(alName.get(i)+"\t"+(cov/(geneBases.length*1.0))+"\n");

			}
			BufferedWriter stat=new BufferedWriter(new FileWriter(outDir+"/stats.out")); 
			stat.write("Occurrences: "+occ+"\nCaused nucleotide changes in genes: "+changes);
			stat.close();
			bw.close();
			bwCov.close();
			
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);

		}

	}

	
	public ArrayList<Info> getCoveredAreas(String id,double fractionCov,int covWindow){
		int stepsize=10;
		double[] cov=coverage.get(id);
		if(cov.length>covWindow){
			ArrayList<Info> areas=new ArrayList<Info>();
			ArrayList<Double> countArray=getCountArray(stepsize,cov);
			double windowCount=initWindowCount(countArray,stepsize,covWindow);
			int start=covWindow/stepsize;
			if(windowCount/(covWindow*1.0)>=fractionCov){
				addArea(0,covWindow,areas);
			}
			for(int i=0;i<countArray.size()-start;i++){
				double subtract=countArray.get(i);
				double add=countArray.get(i+start);
				windowCount=windowCount+add-subtract;
				double fraction=windowCount/(covWindow*1.0);
				if(fraction>=fractionCov){
					addArea(i*stepsize,i*stepsize+covWindow,areas);
				}
			}
			//System.out.println(areas.size());
			return areas;
		}else{
			if(getfractionCoverage(id)>=fractionCov){
				ArrayList<Info> wholeArea=new ArrayList<Info>();
				wholeArea.add(new Info(0,cov.length,""));
				return  wholeArea;
			}else{
				return new ArrayList<Info>();
			}
		}
	}
	
	
	/**
	 * Extracts the nucleotide position in the original fasta file, given the mapping of the read.
	 * @param readID
	 * @param subpos
	 * @param orientation
	 * @param length
	 * @return
	 */
	 int getQueryName(String readID,int subpos,char orientation,int length){
		 int readIDPos=Integer.parseInt(readID);
		 subpos=orientation=='+'?subpos:length-(subpos+1);
		 int position=subpos+readIDPos;

		 return position;
		
	}
	
	 String getStrain(String readIDsplit[]){
		StringBuffer strain=new StringBuffer();
		strain.append(readIDsplit[0]);
		for(int i=1;i<readIDsplit.length-1;i++){
			strain.append("_"+readIDsplit[i]);
		}
		return strain.toString();
	}
	
	public HashMap<String,Integer> getSubstitutionInfo(String refID){
		
		return substitutionInfo.get(refID);
	}
	
	public int getMappedSites(int minCov){
		int sum=0;

		Iterator<Entry<String,double[]>> it=coverage.entrySet().iterator();
		while(it.hasNext()){
			double[] cov=it.next().getValue();
			for(int i=0;i<cov.length;i++){
				if(cov[i]>=minCov){
					sum++;
				}
			}
		}
		return sum;
	}
	

	 void setSubstitution(String readSequence,int mismatchRead,int mismatchRef,String readID,String qualityString,char orientation,int posorig,int flank,int quality,String geneId,boolean subInfo,double weight){
		int lengthSeq=substitutions.get(geneId).length;
		int sub = posorig+mismatchRef-flank;
		if(sub>0&& sub<lengthSeq){
			if((int)(qualityString.charAt(mismatchRead))-33>=quality){
				String base = orientation=='+'?readSequence.charAt(mismatchRead)+"F":DNAmanipulations.reverse(readSequence.charAt(mismatchRead)+"")+"R";
				if (lengthSeq > sub){
					substitutions.get(geneId)[sub]+=weight;
					if(!subInfo){
						//System.out.println(sub+" "+base);
						bases.get(geneId).set(sub, base,weight);
					}else{
						int length=readSequence.length();
						if(!bases.get(geneId).isset(sub)){
							int subID=getQueryName(readID, mismatchRead, orientation,length);
							bases.get(geneId).set(sub, base,subID,weight);
						}else{
							bases.get(geneId).set(sub, base,weight);
						}

					}
					
				}
			}
		}

	}
	 
	 
	 //TODO I think there is something wrong here!!! I can't name a gap position if it does not exist in the query.
	 /**
	  * This method is supposed to set a gap at the designated positions in the reference sequence. This position is determined by adding posorig and mismatchref. 
	  * The position of the gap in the read is needed to calculate the correct name of the position.
	  * @param readSequence
	  * @param mismatchRead
	  * @param mismatchRef
	  * @param readID
	  * @param orientation
	  * @param posorig
	  * @param flank
	  * @param fastaId
	  * @param subInfo
	  * @param weight
	  */
	 void setGap(String readSequence,int mismatchRead,int mismatchRef,String readID,char orientation,int posorig,int flank,String fastaId,boolean subInfo,double weight){
		 int lengthSeq=substitutions.get(fastaId).length;
		 int sub = posorig+mismatchRef-flank;
		 if(sub>0&& sub<lengthSeq){
			 String base = orientation=='+'?"-F":"-R";
			 if (lengthSeq > sub){
				 substitutions.get(fastaId)[sub]+=weight;
				 if(!subInfo){
					 bases.get(fastaId).set(sub, base,weight);
				 }else{
					 int length=readSequence.length();
					 if(!bases.get(fastaId).isset(sub)){
						 int subID=getQueryName(readID, mismatchRead, orientation,length);
						 bases.get(fastaId).set(sub, base,subID,weight);
					 }else{
						 bases.get(fastaId).set(sub, base,weight);
					 }
				 }

			 }
		 }
	 }
	 /**
	  * Sets the coverage/name of a range of bases. Only works if there are no gaps in both read and reference sequence.
	  * 
	  */
	 
	 void setCoverage(int pos,int length,String fastaId,String sequence,String qualityString,int quality,char orientation,boolean subInfo,String readID, double weight){
		int lengthSeq=	coverage.get(fastaId).length;
		double[] cov=coverage.get(fastaId);
		double[] covPos=coveragePos.get(fastaId);
		double[] covNeg=coverageNeg.get(fastaId);
		//System.out.println(readID+" "+pos+" "+length);
		 for (int i = pos; i < pos + length && i<lengthSeq; i++) {
			 int subpos=i-pos;
			 if (i>0&&lengthSeq > i &&qualityString.charAt(subpos)-33>=quality){
				 cov[i]+=weight;

				 if(orientation=='+'){
					 covPos[i]+=weight;
				 }else{
					 covNeg[i]+=weight;
				 }
			 }
			 if(subInfo){
				 
				 addBaseNames(fastaId, i, readID, subpos, orientation, length);

			 }
		 }
		 
	 }
	 /**
	  * Sets the coverage for a single genome position and if gap is set to false then it also sets the name for that position.
	  * 
	  * @param pos
	  * @param readPos
	  * @param fastaId
	  * @param sequence
	  * @param qualityString
	  * @param quality
	  * @param orientation
	  * @param subInfo
	  * @param readID
	  * @param weight
	  * @param gap
	  */
	 void setCoverageSingle(int pos,int readPos,String fastaId,String sequence,String qualityString,int quality,char orientation,boolean subInfo,String readID, double weight,boolean gap,int flank){
		int lengthSeq=	coverage.get(fastaId).length;
		double[] cov=coverage.get(fastaId);
		double[] covPos=coveragePos.get(fastaId);
		double[] covNeg=coverageNeg.get(fastaId);
		int subpos=readPos;
		int sub=pos-flank;
		if (sub>=0&&lengthSeq > sub &&qualityString.charAt(subpos)-33>=quality){
			cov[sub]+=weight;

			if(orientation=='+'){
				covPos[sub]+=weight;
			}else{
				covNeg[sub]+=weight;
			}
		}
		if(subInfo&&!gap){
			addBaseNames(fastaId, sub, readID, subpos, orientation,qualityString.length());
		}
	 }
		 
	 void addBaseNames(String fastaId,int sub,String readID,int subpos,char orientation,int length){
		 int subID=getQueryName(readID, subpos, orientation,length);

		 if(baseNames.get(fastaId)[sub]==null){
				HashMap<Integer,Integer> temp=new HashMap<Integer,Integer>();
				temp.put(subID,1);
				baseNames.get(fastaId)[sub]=temp;
			}else{
				if(baseNames.get(fastaId)[sub].containsKey(subID)){
					int item=baseNames.get(fastaId)[sub].get(subID);
					item++;
					baseNames.get(fastaId)[sub].put(subID,item);
				}
			}
	 }
	
}
