package buildTreesFromSequenceData;


import java.io.*;
import java.util.*;
import java.util.Map.Entry;





import solexa.*;
import solexa.Arrays;
import util.*;

	//input is a folder that contains all the alignment files (strain number with extension sop) and a core gene set
	//alignment files are analyzed with point substitutions one after the other, for each gene polymorphism sites are stored
	//based on coverage core gene set is reduced
	//output is a multi fasta sequence consisting of DNA codons
	public class GetPolymorphismWithGaps extends GetPolymorphismsClass {
		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;




		public GetPolymorphismWithGaps(ArrayList<File> alignmentFiles,ArrayList<String> references,File ReferenceFas,int Flank,int Quality,double PolymorphismThreshold,double fractionCovThreshold,int PerBaseCoverage,boolean subInfo,boolean NoGenes,double gapThreshold,boolean printInvariant,File outFolder,String reference) {
			super(alignmentFiles,references,ReferenceFas,Flank,Quality,PolymorphismThreshold,fractionCovThreshold,PerBaseCoverage, subInfo,NoGenes,printInvariant, outFolder,true,gapThreshold,reference);
			
			

		}
		
		private void checkGapThreshold(double gapThreshold){
			int strainNumber=strains.size();
			int genePolySize=genePoly.size();
			for(int i=0;i<genePolySize;i++){
				Iterator<Entry<Integer,Integer>> it=genePoly.get(i).gaps.entrySet().iterator();
				String geneID=genePoly.get(i).geneID;
				while(it.hasNext()){
					Entry<Integer,Integer> e=it.next();
					int gaps=e.getValue();
					double gapProportion=(1.0*gaps)/strainNumber;
					
					if(gapThreshold<=gapProportion){
						genePoly.get(i).polymorphisms.remove(e.getKey());
						
						//also delete the position in all other strains to save time and space
						Iterator<Entry<String,HashMap<String,HashMap<Integer,Polymorph>>>> it2=strains.entrySet().iterator();
						while(it2.hasNext()){
							Entry<String,HashMap<String,HashMap<Integer,Polymorph>>> e2=it2.next();
							e2.getValue().get(geneID).remove(e.getKey());
						}
					
					}
				}
			}
			
		}
		

		

		


		

		

		
		public void checkCoverageGP(PointSubstitutions pss,String strain){

				HashMap<String,HashMap<Integer,Polymorph>> strainPoly=strains.get(strain);
				int genePolySize=genePoly.size();
				for(int i=0;i<genePolySize;i++){
					String id=genePoly.get(i).geneID;

					ArrayList<Integer> gaps=deleteUnsureRegions(pss, i);

					HashMap<Integer,Polymorph> geneHM=null;
					if(!strainPoly.containsKey(id)){
						geneHM=new HashMap<Integer, Polymorph>();
						strainPoly.put(id,geneHM);
					}else{
						geneHM=strainPoly.get(id);
					}
					int gapsSize=gaps.size();
					for(int j=0;j<gapsSize;j++){
						int pos=gaps.get(j);
						Genes genes=genePoly.get(i);

						geneHM.put(pos,new Polymorph('-'));

						if(genes.gaps.containsKey(pos)){
							int number=genes.gaps.get(pos);
							genes.gaps.put(pos,number+1);
						}else{
							genes.gaps.put(pos,1);
						}

					}
				}
				checkGapThreshold(gapThreshold);

		}
		

		
		public File writeSequences(){
			File fas=new File(outFolder+"/polymorphisms.fas");
			try{
				BufferedWriter bw=new BufferedWriter(new FileWriter(fas));
				Iterator<Entry<String,HashMap<String,HashMap<Integer,Polymorph>>>> it=strains.entrySet().iterator();

				while(it.hasNext()){

					Entry<String,HashMap<String,HashMap<Integer,Polymorph>>> e=it.next();
					String strainIDIntern=e.getKey();
					String strainIDExtern=getStrain(strainIDIntern);
					BufferedWriter bwDet=new BufferedWriter(new FileWriter(new File(outFolder+"/"+strainIDExtern+"details.txt")));
					BufferedWriter bwGenePoly=new BufferedWriter(new FileWriter(new File(outFolder+"/"+strainIDExtern+"_genePoly.txt")));
					bwGenePoly.write("Strain\tContig\tGene\tRefPos\tSyn\tNonSyn\tratio\tgeneLength\tpoly/geneLength\n");
					bwDet.write("Strain\tContig\tGene\tRefGeneStartPos\tRefGenePos\tRefGenomePos\tOrig\tPoly\n");
					bw.write(">"+strainIDExtern+"\n");
					HashMap<String,HashMap<Integer,Polymorph>> strainPoly=e.getValue();
					int posInAlignment=1;
					for(int i=0;i<genePoly.size();i++){
						HashMap<Integer,Polymorph> polies=strainPoly.get(genePoly.get(i).geneID);
						Iterator<Entry<Integer,Character>> it2= genePoly.get(i).polymorphisms.entrySet().iterator();
						int syn=0;
						int nonsyn=0;
						String contig=genePoly.get(i).contig;
						
						while(it2.hasNext()){
							if(polies==null){
								System.err.println(genePoly.get(i).geneID+" "+strainIDExtern);
							}
							Entry<Integer,Character> e2=it2.next();
							int pos=e2.getKey();
							char origBase=e2.getValue();
							if(polies.containsKey(pos)){
								int GenomePos=genePoly.get(i).orient=='+'?pos+genePoly.get(i).pos-1:(genePoly.get(i).length-pos)+genePoly.get(i).pos;
								bwDet.write(strainIDExtern+"\t"+contig+"\t"+genePoly.get(i).geneID+"\t"+genePoly.get(i).pos+"\t"+pos+"\t"+GenomePos+"\t"+origBase+"\t"+polies.get(pos).base+"\t"+posInAlignment);
								if(genes){
									int frame=(pos-1)%3;
									String origCodon=""+genePoly.get(i).polymorphisms.get(pos-frame)+""+genePoly.get(i).polymorphisms.get(pos-frame+1)+""+genePoly.get(i).polymorphisms.get(pos-frame+2);
									origCodon+="="+DNAmanipulations.translate(origCodon, DNAmanipulations.code());
									char one=polies.containsKey(pos-frame)?polies.get(pos-frame).base:genePoly.get(i).polymorphisms.get(pos-frame);
									char two=polies.containsKey(pos-frame+1)?polies.get(pos-frame+1).base:genePoly.get(i).polymorphisms.get(pos-frame+1);
									char three=polies.containsKey(pos-frame+2)?polies.get(pos-frame+2).base:genePoly.get(i).polymorphisms.get(pos-frame+2);
									//ADD NAME INFO
									String newCodon=""+one+""+two+""+three;
									newCodon+="="+DNAmanipulations.translate(newCodon, DNAmanipulations.code());
									bwDet.write("\t"+origCodon+"\t"+newCodon+"\n");
									if(DNAmanipulations.translate(origCodon, DNAmanipulations.code()).equals(DNAmanipulations.translate(newCodon, DNAmanipulations.code()))){
										syn++;
									}else{
										nonsyn++;
									}
								}
								else bwDet.write("\n");
								posInAlignment++;
								bw.write(polies.get(pos).base);
							}else {
								posInAlignment++;
								bw.write(origBase);
							}
						}
						if(genes)bwGenePoly.write(strainIDExtern+"\t"+contig+"\t"+genePoly.get(i).geneID+"\t"+genePoly.get(i).pos+"\t"+syn+"\t"+nonsyn+"\t"+(syn/(1.0*nonsyn))+"\t"+genePoly.get(i).length+"\t"+((syn+nonsyn+1.0)/genePoly.get(i).length)+"\n");
					}
					bwDet.close();
					bwGenePoly.close();
					bw.write("\n");
				}

				bw.close();
			}catch(IOException e){
				e.printStackTrace();
			}
			return fas;
		}
		

		/**
		 * Add SNP sites (invar = false) or SNP and conserved sites (invar = true) to the overall site collection as well as strain SNP collection if a set of conditions are fulfilled.
		 */
			

		    public  HashMap<String /*geneID*/,HashMap<Integer/*position*/,Polymorph/*poly*/>> adjustGP(PointSubstitutions pss,boolean subInfo,String name){
		    	HashMap<String /*geneID*/,HashMap<Integer/*position*/,Polymorph/*poly*/>> strainPolymorphism=new HashMap<String /*geneID*/,HashMap<Integer/*position*/,Polymorph/*poly*/>>();
		    	try{
		    		BufferedWriter bw=new BufferedWriter(new FileWriter(outFolder+"/"+name+"deletedSites.txt"));
		    		int genePolySize=genePoly.size();

		    		//go through all genes
		    		for(int i=0;i<genePolySize;i++){
		    			String id=genePoly.get(i).geneID;
		    			Arrays geneBases=pss.getBases(id);
		    			String gene=pss.getGene(id);
		    			double[] geneCoverage=pss.getCoverage(id); 
		    			if(!strainPolymorphism.containsKey(id)){			
		    				strainPolymorphism.put(id, new HashMap<Integer, Polymorph>());
		    			}
		    			//go through all sites
		    			for(int j=1;j<geneCoverage.length;j++){

		    				//if all conditions for a SNP or site are fulfilled 
		    				//the base coverage has to be greater than the min coverage
		    				if(geneCoverage[j]>=perBaseCoverage){
		    					char base=geneBases.getMaxNucForward(j);
		    					double numMajorPoly=geneBases.numMajorPolymorphism(j);
		    					double ratioPoly=numMajorPoly/(geneCoverage[j]*1.0);

		    					//	 SNP ratio > min. polymorphism threshold/ if invar then it is also ok if the SNP ration is smaller than 1-polyT (few errors) 
		    					if ((ratioPoly >= polymorphismThreshold ||(printInvariant&& ratioPoly<= (1-polymorphismThreshold))) ){
		    						//and the site is not already stored 
		    						if(!genePoly.get(i).polymorphisms.containsKey(j)){
		    							int frame=genes?(j-1)%3:0;
		    

		    							//add it to the overall polymorphism collection
		    							genePoly.get(i).polymorphisms.put(j-frame,gene.charAt(j-1+flank-frame));
		    							if(genes){
		    								genePoly.get(i).polymorphisms.put(j-frame+1,gene.charAt(j-1+flank-frame+1));

		    								genePoly.get(i).polymorphisms.put(j-frame+2,gene.charAt(j-1+flank-frame+2));
		    							}
		    						}
		    						//as well as to the strain SNP collection (but only if it is actually a SNP!!!)
		    						if(numMajorPoly/(geneCoverage[j]*1.0) >= polymorphismThreshold){

		    							if(!subInfo)strainPolymorphism.get(id).put(j, new Polymorph(base));

		    							//if we intend to merge SNP files from multiple references then we also need to store information on where the SNP came from
		    							else strainPolymorphism.get(id).put(j, new Polymorph(base,getMajority(geneBases.getQueryID(j))));
		    						}
		    					}else {
		        					bw.write(id+"\t"+j+"\t"+base+"\t"+geneCoverage[j]+"\t"+numMajorPoly+"\t"+ratioPoly+"\n");
		        				}
		    				}
		    			}
		    		}
		    		bw.close();
				}catch(IOException e){
					e.printStackTrace();
					System.exit(-1);
				}
				return strainPolymorphism;
			}

	    
/*	     public HashMap<String geneID,HashMap<Integerposition,Polymorphpoly>> adjustGP(PointSubstitutions pss,boolean subInfo,String name){
			HashMap<String geneID,HashMap<Integerposition,Polymorphpoly>> strainPolymorphism=new HashMap<String geneID,HashMap<Integerposition,Polymorphpoly>>();
			try{
				BufferedWriter bw=new BufferedWriter(new FileWriter(outFolder+"/"+name+"deletedSites.txt"));

				int genePolySize=genePoly.size();
				for(int i=0;i<genePolySize;i++){
					String id=genePoly.get(i).geneID;
					Arrays geneBases=pss.getBases(id);
					String gene=pss.getGene(id);
					double[] geneCoverage=pss.getCoverage(id);
					ArrayList<Info> wholeArea=new ArrayList<Info>();
					wholeArea.add(new Info(0,geneCoverage.length,""));
					//ArrayList<Info> coveredAreas=covWindow>0?pss.getCoveredAreas(id,fractionCovThreshold,covWindow):wholeArea;
					int coverSize=wholeArea.size();
					for(int k=0;k<coverSize;k++){
						if(!strainPolymorphism.containsKey(id)){			
							strainPolymorphism.put(id, new HashMap<Integer, Polymorph>());
						}
						TreeMap<Integer,Character> genePolymorphismPositions=genePoly.get(i).polymorphisms;
						Info area=wholeArea.get(k);
						int start=area.getStart();
						int end=area.getEnd();
						for(int j=start;j<end;j++){

							double numMajorPoly=geneBases.numMajorPolymorphism(j);
							if(geneCoverage[j]>=perBaseCoverage){
								char base=geneBases.getMaxNucForward(j);
								double ratioPoly=numMajorPoly/(geneCoverage[j]*1.0);

								if ( ratioPoly>= polymorphismThreshold||(printInvariant&&ratioPoly <= (1-polymorphismThreshold))) {
									if(!genePolymorphismPositions.containsKey(j)){
										int frame=genes?(j-1)%3:0;
										genePolymorphismPositions.put(j-frame,gene.charAt(j-1+flank-frame));
										if(genes){
											genePolymorphismPositions.put(j-frame+1,gene.charAt(j-1+flank-frame+1));

											genePolymorphismPositions.put(j-frame+2,gene.charAt(j-1+flank-frame+2));
										}
									}
									if(ratioPoly >= polymorphismThreshold){
										if(!subInfo)strainPolymorphism.get(id).put(j, new Polymorph(base));
										else strainPolymorphism.get(id).put(j, new Polymorph(base,geneBases.getQueryID(j)));
									}
								}else {
									bw.write(id+"\t"+j+"\t"+base+"\t"+geneCoverage[j]+"\t"+numMajorPoly+"\t"+ratioPoly+"\n");
								}
							}
						}
					}
				}
				bw.close();
			}catch(IOException e){
				e.printStackTrace();
				System.exit(-1);
			}
			return strainPolymorphism;
		}
*/
		
	}


