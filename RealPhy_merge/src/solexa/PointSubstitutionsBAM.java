package solexa;

import java.io.*;
import java.util.*;


import net.sf.samtools.*;



public class PointSubstitutionsBAM extends PointSubstitutions {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	public static void main(String args[]){
		File refseq=new File("/home/frederic/Basel/randomAlignments/test/test2/sequences/S21.fas");
		PointSubstitutionsBAM pss=new PointSubstitutionsBAM(refseq, 0, new File("/home/frederic/Basel/randomAlignments/test/test2/S21/alignOut_NoGenes/S12_50fasta.bam"), 20, false);
		pss.writeArray(0.0, refseq, new File("/home/frederic/Basel/randomAlignments/test/test2/S21/alignOut_NoGenes/"), 0);
		
	}
	
	public void setSubstitution(String MD,String cigar,String sequence,String readID,String qualityString,char orientation,int posorig,int flank,int quality,String fastaId,boolean subInfo,double weight){

		ArrayList<Integer> posMismatchesRef=getMismatchesRef(MD);
		ArrayList<Integer> posMismatchesRead;
		if(cigar.contains("I")){
			ArrayList<Integer> refDeletions=getRefDeletions(cigar);
			posMismatchesRead=getMismatchesRead(MD,refDeletions);
		}else if(cigar.contains("D")){
			posMismatchesRead=getMismatchesRead(MD);

		}else{
			posMismatchesRead=posMismatchesRef;
		}
		
		
		setSubstitutions( posMismatchesRead, posMismatchesRef, sequence, readID, qualityString, orientation, posorig, flank, quality, fastaId, subInfo,weight);
	}
	
	public void setGaps(String cigar,String sequence,String readID,char orientation,int posorig,int flank,String fastaId,boolean subinfo,double weight){
		ArrayList<Integer> gapsRead=getGapsRead(cigar);
		ArrayList<Integer> gapsRef=getGapsRef(cigar,orientation);
		setGaps(gapsRead,gapsRef,sequence,readID,orientation,posorig,flank,fastaId,subinfo,weight);
	}
	
	private ArrayList<Integer> getGapsRead(String cigar){
		String pos[]=cigar.split("[^0-9]+");
		String edit[]=cigar.split("[0-9]+");
		ArrayList<Integer> gapsPos=new ArrayList<Integer>();
		int start=0;
		for(int i=1;i<pos.length;i++){
			int count=Integer.parseInt(pos[i-1]);

			if(edit[i].charAt(0)=='D'){
				for(int j=0;j<count;j++){
					gapsPos.add(start);
				}
			}else if(edit[i].charAt(0)=='M'||edit[i].charAt(0)=='I'){

				start=start+count;
				
			}
		}
		return gapsPos;
	}	
	private ArrayList<Integer> getGapsRef(String cigar,char orient){
		String pos[]=cigar.split("[^0-9]+");
		String edit[]=cigar.split("[0-9]+");
		ArrayList<Integer> gapsPos=new ArrayList<Integer>();
		int start=0;
		for(int i=1;i<pos.length;i++){
			int count=Integer.parseInt(pos[i-1]);

			if(edit[i].charAt(0)=='D'){
				for(int j=0;j<count;j++){
					gapsPos.add(start+j);
				}
				start=start+count;
			}else if(edit[i].charAt(0)=='M'){

				start=start+count;
				
			}
		}
		return gapsPos;
	}
	
	public void setCoverageCigar(StringBuffer regionsRef,StringBuffer regionsRead,int length,int pos,String fastaId,String sequence,String qualityString,int quality,char orientation,boolean subInfo,String readID,double weight,int flank){

		int size=regionsRef.length();
		int readPos=-1;
		int posorig=pos;
		boolean gap=true;
		if(orientation=='+'){
			for(int i=0;i<size;i++){
				if(regionsRead.charAt(i)=='m'){
					readPos++;
					gap=false;
				}
				if(regionsRef.charAt(i)=='m'){
					setCoverageSingle(pos, readPos, fastaId, sequence,qualityString, quality, orientation, subInfo, readID,weight,gap,flank);
					if(gap){
						setGap(sequence, readPos, pos-posorig, readID, orientation, posorig, flank, fastaId, subInfo, weight);
					}
					pos++;
				}

				gap=true;
			} 
		}else{
			for(int i=size-1;i>=0;i--){
				if(regionsRead.charAt(i)=='m'){
					readPos++;
					gap=false;
				}
				if(regionsRef.charAt(i)=='m'){
					setCoverageSingle(pos, readPos, fastaId, sequence,qualityString, quality, orientation, subInfo, readID,weight,gap,flank);
					if(gap){
						setGap(sequence, readPos, pos-posorig, readID, orientation, posorig, flank, fastaId, subInfo, weight);
					}
					pos++;
					
				}

				gap=true;
			} 
		}

	}
	
	public void setSubstitutions(ArrayList<Integer> posMismatchesRead,ArrayList<Integer> posMismatchesRef,String sequence,String readID,String qualityString,char orientation,int posorig,int flank,int quality,String fastaId,boolean subInfo, double weight){
		for(int j=0;j<posMismatchesRead.size();j++){
			int mismatchRead=posMismatchesRead.get(j)-1;
			int mismatchRef=posMismatchesRef.get(j)-1;
			
			setSubstitution(sequence, mismatchRead,mismatchRef, readID, qualityString, orientation, posorig, flank, quality, fastaId,  subInfo,weight);

		}
	}

	public void setGaps(ArrayList<Integer> posGapsRead,ArrayList<Integer> posGapsRef,String sequence,String readID,char orientation,int posorig,int flank,String fastaId,boolean subInfo, double weight){
		for(int j=0;j<posGapsRef.size();j++){
			int gapRead=posGapsRead.get(j);
			int gapRef=posGapsRef.get(j);

			setGap(sequence,gapRead,gapRef, readID, orientation, posorig, flank, fastaId,  subInfo,weight);

		}
	}

		
		public PointSubstitutionsBAM(File RefSeq,int flank,File AlignmentFile,int quality,int fold,boolean subInfo){
			super( RefSeq, flank, AlignmentFile, quality, fold, subInfo);
		}
		public PointSubstitutionsBAM(File RefSeq,int flank,File AlignmentFile,int quality,boolean subInfo){
			this(RefSeq,flank,AlignmentFile,quality,1,subInfo);
		}
		@Override
		void read(int quality, int flank, int fold) {
				SAMFileReader sfr=new SAMFileReader(alignmentFile,  true);

				SAMRecordIterator sri=sfr.iterator();
				
				
				double p=1.0/fold;
				SAMRecord srnext=null;
				while (sri.hasNext()) {
					if(p<1&&Math.random()>p){
						sri.next();	
						continue;
					}
					ArrayList<SAMRecord> sams=new ArrayList<SAMRecord>();
					sams.add(srnext!=null?srnext:(srnext=sri.next()));
					SAMRecord sr=srnext;
					int AS=(Integer)sr.getAttribute("AS");
					while(sri.hasNext()&&sr.getReadName().equals((srnext=sri.next()).getReadName())){//&&srnext.getFlags()>=256){
						int ASnext=(Integer)srnext.getAttribute("AS");
						boolean readPair=sr.getReadPairedFlag()&&sr.getFirstOfPairFlag()&&srnext.getSecondOfPairFlag();
						if(AS==ASnext||readPair)sams.add(srnext);
						sr=srnext;
					}
					
					analyseAll(sams,quality,flank);
					
						
					
				}
				sri.close();
				sfr.close();
		
			
		}
		
		private void analyseAll(ArrayList<SAMRecord> sams,int quality,int flank){
			double weight=1/(sams.size()*1.0);
			//System.out.println(sams.size()+" "+weight);
			for(int i=0;i<sams.size();i++){
				
				analyseLine(sams.get(i), flank, quality, weight);
			}
		}
		
		private void analyseLine(SAMRecord sr,int flank,int quality,double weight){
			if(sr.getReadPairedFlag() && sr.getProperPairFlag()){
				weight=weight*2;
			}
			if(weight>1){
				System.err.println("Weight is too large (>1): "+weight+" for read "+sr.getReadName()+" in file "+this.alignmentFile+".");
				System.exit(-1);
			}
			String cigar=sr.getCigarString();
			int length=sr.getReadLength();
			String qualityString=sr.getBaseQualityString();
			String fastaId=sr.getReferenceName();
			String sequence=sr.getReadString();
			char orientation=sr.getReadNegativeStrandFlag()?'-':'+';
			String readID=sr.getReadName();
			int posorig=sr.getAlignmentStart();;
			//System.out.println(posorig);
			int pos=posorig;

			if(pos<=flank){
				length=(pos+length)-flank-1;
				pos=1;
			}else{
				pos=pos-flank;
			}
			Object XM=sr.getAttribute("XM");
			//System.out.println(cigar+" "+ length+" "+pos+" "+fastaId+" "+sequence+" "+qualityString+" "+quality+" "+orientation+" "+subInfo+" "+readID);
			
			
			StringBuffer regionsRef=new StringBuffer();
			StringBuffer regionsRead=new StringBuffer();
			
			getMatchRegionsCigar(cigar,regionsRead,regionsRef);
			//System.out.println(sr);
			if(!cigar.equals(length+"M")){
				setCoverageCigar(regionsRef,regionsRead, length, pos, fastaId, sequence, qualityString, quality, orientation, subInfo, readID,weight,flank);
			}else{
				setCoverage(pos, length, fastaId, sequence,qualityString, quality, orientation, subInfo, readID,weight);

			}
			//setGaps(cigar,sequence,readID,orientation,posorig,flank,fastaId,subInfo,weight);
			if(XM!=null&&(Integer)XM>0){
				Object MD=sr.getAttribute("MD");

				//System.out.println(MD);
				//System.out.println(MD.toString()+" "+ cigar+" "+ sequence+" "+ readID+" "+ qualityString+" "+ orientation+" "+ posorig+" "+ flank+" "+ quality+" "+ fastaId+" "+ subInfo);
				setSubstitution(MD.toString(), cigar, sequence, readID, qualityString, orientation, posorig, flank, quality, fastaId, subInfo,weight);
			}
		}

		
	}
