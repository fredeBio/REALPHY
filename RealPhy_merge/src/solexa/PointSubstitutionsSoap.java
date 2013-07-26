package solexa;
import java.io.*;


public class PointSubstitutionsSoap extends PointSubstitutions implements Serializable{
	/**
	 * Deprecated!!!
	 */
	private static final long serialVersionUID = 1L;
	public PointSubstitutionsSoap(File RefSeq,int flank,File AlignmentFile,int quality,int fold,boolean subInfo){
		super( RefSeq, flank, AlignmentFile, quality, fold, subInfo);
	}
	public PointSubstitutionsSoap(File RefSeq,int flank,File AlignmentFile,int quality,boolean subInfo){
		this(RefSeq,flank,AlignmentFile,quality,1,subInfo);
	}
	public static void main(String[] args) {

		File soap=new File(args[0]);
		File RefSeq = new File(args[1]);
		double threshold = Double.parseDouble(args[2]);
		File outDir = new File(args[3]);
		int quality=Integer.parseInt(args[4]);
		int flank=Integer.parseInt(args[5]);

		PointSubstitutionsSoap PS=new PointSubstitutionsSoap(RefSeq, flank, soap, quality,false);
		

		PS.writeArray( threshold,  RefSeq, outDir,flank);

	}
	


	

	
	 void read( int quality,int flank,int fold) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(alignmentFile));
			String line = "";
			double p=1.0/fold;
			while ((line = br.readLine()) != null) {
				if(Math.random()>p)continue;
				String data[] = line.split("\\s+");
				if(Integer.parseInt(data[3])>1)continue; 
				String fastaId=data[7];
				char orientation=data[6].charAt(0);
				int posorig = Integer.parseInt(data[8]);
				int length = Integer.parseInt(data[5]);
				int pos=posorig;
				String sequence = data[1];
				String readID=data[0];
				String qualityString=data[2];
				if(pos<=flank){
					length=(pos+length)-flank-1;
					pos=1;
				}else{
					pos=pos-flank;
				}
				if(!coverage.containsKey(fastaId)){
					System.err.println("The reference name in the soap alignment file is "+fastaId+". There is no such reference name in the genbank or fasta file.\nExiting!");
					//System.err.println(coverage.keySet().toArray().length);
					System.exit(-1);
				}
				
				setCoverage(pos, length, fastaId, sequence,qualityString, quality, orientation, subInfo, readID,1.0);
				int numSubs=Integer.parseInt(data[9]);
				//
				//Exclude non-unique hits/ Exclude repeats
				//
				for(int j=0;j<numSubs;j++){
					//setSubstitution(data, posorig, flank, quality, fastaId, j,subInfo);
					 String[] split=data[10+j].split("[^0-9]+");

					 int readSubPos=Integer.parseInt(split[1]);

					 setSubstitution(sequence, readSubPos,readSubPos, readID, qualityString, orientation, posorig, flank, quality, fastaId,  subInfo,1.0);

				}
				//System.exit(1);

			}
			br.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}

	}
	




	

		
/*	public void readSeq(File ref, int flank) {
		try {
			BufferedReader is = new BufferedReader(new FileReader(ref));
			String line = "";
			StringBuilder sb = new StringBuilder();
			String name = "";
			while ((line = is.readLine()) != null) {
				if (!line.startsWith(">"))
					sb.append(line.trim());
				else {
					if (name.length() > 0) {
						alName.add(name.split("\\s+")[0]);
						alLength.add(sb.length()-2*flank);
						sb = new StringBuilder();
					}
					name = line.substring(1);

				}
			}
			alName.add(name.split("\\s+")[0]);
			alLength.add(sb.length());
			is.close();
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}*/
//	private ArrayList<String> getAlt(String triplet,int pos, int frame,Arrays bases,boolean complement,int minCov){
//		ArrayList<String> alts=new ArrayList<String>();
//		char orient=complement?'R':'F';
//		if(bases.get(pos, "A"+orient)>minCov){
//			StringBuilder temp=new StringBuilder(triplet);
//			temp.setCharAt(frame, 'A');
//			alts.add(temp.toString());
//		}
//		if(bases.get(pos, "T"+orient)>minCov){
//			StringBuilder temp=new StringBuilder(triplet);
//			temp.setCharAt(frame, 'T');
//			alts.add(temp.toString());
//		}
//		if(bases.get(pos, "C"+orient)>minCov){
//			StringBuilder temp=new StringBuilder(triplet);
//			temp.setCharAt(frame, 'C');
//			alts.add(temp.toString());
//		}
//		if(bases.get(pos, "G"+orient)>minCov){
//			StringBuilder temp=new StringBuilder(triplet);
//			temp.setCharAt(frame, 'G');
//			alts.add(temp.toString());
//		}
//		return alts;
//	}

}
