package buildTreesFromSequenceData;

import java.io.*;
import java.util.*;




import rcaller.*;
import util.phylogenetics.RunTreePrograms;
import util.*;

/*
REALPHY, a program to automatically reconstruct phylogenetic trees from short sequence reads. 

Copyright (C) 2014 Frederic Bertels
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
public class RealPhy {
	
	static String version="v1.12";
	
	public static  final String fasExt[]=new String[]{"fas","fa","fasta","fna","fsa"};
	public static final String gbkExt[]=new String[]{"gbk","gb"};
	public static final String fastqExt[]=new String[]{"fastq","fq"};
	public static final String gzExt[]=new String[]{"fastq.gz"};

	File sequenceFolder=null;
	File outFolder=null;
	File alignmentFolder=null;
	File masterOutFolder=null;

	File mergeFolder=null;
	ArrayList<String> reflist=new ArrayList<String>();
	boolean clean;
	File config=null;
	int seedLength=20;
	HashMap<String,Object> arguments=new HashMap<String,Object>();
	HashMap<String,String> path=new HashMap<String, String>();
	boolean references=false;
	boolean delete=false;
	public final int bowtie2=1;
	public final int soap2=0;
	int aligner=bowtie2;
	int perBaseCov=10;
	boolean runAlignment=false;
	File polymorphismsOutFolder;
	int flank=0;
	File cutFolder;
	Boolean gaps=false;
	Double gapThreshold=0.0;
	File treeOptions=null;
	File bowtieOptions=null;

	ArrayList<String> refs;
	public static void main(String args[])throws RealphyException{
		long l=System.currentTimeMillis();
		if(args.length>1){
			RealPhy bt=new RealPhy(args);
			if((Boolean)bt.arguments.get("merge")){
				
				bt.mergeAnalyses();
			
			}
			else {
				bt.runAnalyses();
			}
			bt.analyseDeletedSites();
		}else{
			if(args.length>0&&args[0].equals("-version")){
				printVersion();
			}else{
			printHelp();
			}
		}
		//Cut sequences
		long l2=System.currentTimeMillis();
		double minutes=(l2-l)/60000.0;
		System.out.println("Program execution took "+minutes+" minutes.");
	}
	
	public void analyseDeletedSites(){
		String RscriptExecutable=getRscript();
		double polT=0.9;
		if(RscriptExecutable!=null){
			AnalyseDeletedSites ads=new AnalyseDeletedSites(perBaseCov,polT,polymorphismsOutFolder);
			File deletedSitesHist=new File(polymorphismsOutFolder+"/deletedSitesHist.jpg");
			File deletedSitesAbovePol=new File(polymorphismsOutFolder+"/deletedSitesAbovePolT"+polT+".txt");
			ads.plotHistograms(deletedSitesHist, new File(RscriptExecutable));
			ads.writeDeletedAbovePolT(deletedSitesAbovePol);
		}
	}
	public static void printVersion(){
		System.out.println("This is REALPHY version "+version+".");
		System.exit(0);
	}
	public static void printHelp(){
		System.out.print("Usage:\n" +
				"java -Xmx[available RAM in MB]m -jar RealPhy_"+version+".jar [Sequence folder] [Output folder] [Options]\n" +
				"Sequence folder needs to contain fasta files ending with .fas, .fna, .fasta or .fa, genbank files ending in .gbk or .gb and short read files in fastq format ending in .fastq or fastq.gz.\n" +
				"The output folder needs to contain a file called \"config.txt\", which contains information about the location of the required executables such as bowtie2.\n\n"+
				"Options:\n" +
				"-readLength [integer] default=50 Possible values: Integer greater than 20; Size of fragments that are to be produced from the FASTA/GBK input sequences. With longer read lengths the mapping will take longer but will also map more divergent sequences.\n" +
				"-quality [integer] default=20; Possible values: Integer value between 0 and 41 that corresponds to quality values in fastq files (PHRED+33). Bases with values below thresold in fastq file will not be considered (fasta files will be converted into fastq files with a quality of 20).\n" +
				"-polyThreshold [double] default=0.95; Possible values: Double value between 0 and 1.  Polymorphisms that occur at lower frequency than the specified threshold at any given position of the alignment will not be considered.\n"+
				"-perBaseCov [integer] default=10; Possible values: Integer greater than or equal to 10.  Polymorphisms will only be extracted for regions that are covered by more than the set threshold of reads.\n" +
				"-ref [sequence file name (without extension or path!)] default=random; Possible values: The file name of a sequence data set without the extension (.fas, .fasta, .fa, .fna, .fastq, .fastq.gz, .gb or .gbk). Sets the reference sequence.\n" +
				"-root [sequence file name (without extension or path!)] default=random; Possible values: The file name of a sequence data set without the extension (.fas, .fasta, .fa, .fna, .fastq, .fastq.gz, .gb or .gbk).  Specifies the root of the tree.\n" +
				//"-test [boolean] default=false; Possible values: true/false; column comparison tool for multiple different references.\n" +
				"-refN [sequence file name (without extension or path!)] where N is the n-th reference genome; default=not set; Possible values: The file name of a sequence data set without the extension (.fas, .fasta, .fa, .fna, .fastq, .fastq.gz, .gb or .gbk).\n" +
				"-genes If set then genes (CDS) are extracted from a given genbank file.\n" +
				"-gapThreshold [double] default=0; specifies the proportion of input sequences that are allowed to contain gaps in the final sequence alignment (i.e. if set to 0.2 at most 20% of all nucleotides in each final alignment column are allowed to be gaps).\n" +
				"-clean/-c If set then the whole analysis will be rerun and existing data will be overwritten!\n" +
				"-treeBuilder [integer] default=4;\n" +
				"   0=Do not build a tree;\n" +
				"   1=treepuzzle; \n" +
				"   2=raxml\n" +
				"   3=max. parsimony (dnapars)\n" +
				"   4=PhyML (default)\n" +
				"-quiet/-q If set then it suppresses any program output except for errors or warnings.\n" +
				"-varOnly/-v If set then homologous positions that are conserved in all input sequences are put out. If set then the reconstructed tree may be wrong.\n" +
				//"-aligner [integer] default=1;\n" +
				//"	0=SOAP2\n" +
				//"	1=bowtie2\n" +
				"-seedLength [integer] default=22 Possible values: Integer between 4 and 32; specifies k-mer length in bowtie2.\n" +
				"-suffix [string] default=not set; appends a suffix to the reference output folder.\n" +
				"-d/-delete If this option is set then all alignment output files and sequence cut files will be deleted after the program is terminated.\n" +
				"-merge/-m If this option is set multiple references are selected, a final polymorphism file will be generated which combines all polymorphism files for all references. \n" +
				"-gaps/-g If this option is set. The gapThreshold is automatically set to 100%, unless a different gapThreshold is specified.\n" +
				"-config [string] this specifies the location of the config.txt. If not set it is assumed that the config.txt is in the working directory.\n" +
				"-treeOptions [text file] This option allows the user to provide command line parameters to the tree program in the first line of a given text file.\n" +
				"-bowtieOptions [text file] This option allows the user to provide command line parameters to bowtie2 in the first line of a given text file.\n" +
				"-h/-help Shows this help.\n" +
				"-version Prints the program's version.\n" );
				
		//" -evaluateOverlaps default=false; Possible values: true/false; For test purposes. Tests whether polymorphisms that are only found for one of two compared reference genomes are really the result of one sequence containing more than two polymorphisms.\n" 
				//"-cutGenes default=false; boolean true/false");
	}
	private static final Map<String, Boolean> novalueArgs;
    static
    {
        novalueArgs = new HashMap<String, Boolean>();
        novalueArgs.put("h", true);
        novalueArgs.put("help", true);
        
        novalueArgs.put("gaps", true);
        novalueArgs.put("g", true);

        novalueArgs.put("merge", true);
        novalueArgs.put("m", true);

        novalueArgs.put("d", true);
        novalueArgs.put("delete", true);

        novalueArgs.put("varOnly", true);
        novalueArgs.put("v", true);

        novalueArgs.put("quiet", true);
        novalueArgs.put("q", true);

        novalueArgs.put("clean", true);
        novalueArgs.put("c", true);

        novalueArgs.put("genes", true);

        novalueArgs.put("version", true);


    }
    
	private String getRscript(){
		String RscriptExecutable=path.get("Rscript");
		if(RscriptExecutable==null){
			System.err.println("Cannot find Rscript executable. Please set \"Rscript\" in config.txt.");
			return null;
		}else{
			return RscriptExecutable;
		}
		
	}
	public RealPhy(String[] args){
		initDefault();
		setVariables(args);
		checkQuiet();
		checkSequenceFiles();
		aligner=(Integer)arguments.get("aligner");
		seedLength=(Integer)arguments.get("seedLength");
		flank=(Boolean)arguments.get("genes")==false?0:(Integer)arguments.get("readLength");
		clean=(Boolean)arguments.get("clean");
	}
	
	private void checkQuiet(){
		if((Boolean)arguments.get("quiet")){
			System.setOut(new PrintStream(new OutputStream(){
				public void write(int b) {
				}
			}));
		}
	}
	
	static String getId(String name){
		String[] split=name.split("\\.");
		StringBuffer sb=new StringBuffer();
		
		sb.append(split[0]);
		int num=hasExtension(name, gzExt)?2:1;
		for(int i=1;i<split.length-num;i++){
			sb.append("."+split[i]);
		}
		return sb.toString();
	}
	
	private void checkSequenceFiles(){
		File[] list=sequenceFolder.listFiles();
		HashMap< String,Boolean> hm=new HashMap<String, Boolean>();
		//HashMap< String,Boolean> hm10=new HashMap<String, Boolean>();
		for(int i=0;i<list.length;i++){
			File file=list[i];
			String name=file.getName();
			if(hasExtension(name,gzExt)||hasExtension(name,fastqExt)||hasExtension(name, fasExt)||hasExtension(name, gbkExt)){
				String id=getId(name);
				//String name10=name.length()>10?name.substring(0,10):name;
				if(!hm.containsKey(id)){//&&!hm10.containsKey(name10)){
				//hm10.put(name,true);
					hm.put(id, true);
				}
				else System.err.println("WARNING: "+file+" has the same name or the same name as another file in the sequence folder.\n If the two file names only differ in the extension (.fas, .fasta, .fa, .fna, .fastq, .fastq.gz, .gb or .gbk) only one of the files will be used for further analyses.");
			}
		}
	}
	
	private void clean(){
		
		if(clean){
			
				
			FileHandler.deleteFolder(cutFolder);
			FileHandler.deleteFolder(outFolder);
		}
	}
	
	public void runAnalyses(){
		refs=!references?getReference():getReferences();
		if(refs.size()==0){
			if((Boolean)arguments.get("genes")){
				System.err.println("ERROR: At least one .gbk file is required in the sequence folder: "+sequenceFolder+" (Option -genes is set.)");
				System.exit(-1);
			}else{
				System.err.println("ERROR: At least one genbank or fasta file is required in the sequence folder: "+sequenceFolder);
				System.exit(-1);
			}
		}
		for(int i=0;i<refs.size();i++){
			outFolder=new File(masterOutFolder+"/"+refs.get(i)+arguments.get("suffix"));
			clean();
			if(!outFolder.exists())outFolder.mkdir();
			String genes=(Boolean)arguments.get("genes")==false?"_NoGenes":"";
			String varOnly=(Boolean)arguments.get("varOnly")==false?"":"_noInvar";
			String base=outFolder+"/PolySeqOut"+genes+varOnly+"/";

			polymorphismsOutFolder=new File(base);
			if(!polymorphismsOutFolder.exists()){
				polymorphismsOutFolder.mkdir();
			}

			
			alignmentFolder=!(Boolean)arguments.get("genes")?new File(outFolder+"/alignOut_NoGenes/"):new File(outFolder+"/alignOut/");

			
			System.out.println("Preparing reference...");
			File core=!(Boolean)arguments.get("genes")==true?moveRefToCore(refs.get(i)):determineCore(refs.get(i));
			ArrayList<File> missingAlignmentFiles=getAlignments(false);
			if(missingAlignmentFiles.size()>0){
				System.out.println("Cut fasta sequence files into fragments...");
				ArrayList<File> cuts=cutSequences();
				System.out.println("Running alignment program to align sequence fragments to reference...");
				runAlignmentProgram(core,cuts);
				FileHandler.deleteFolder(cutFolder);
			}
			ArrayList<File> alignmentFiles=getAlignments(true);
			System.out.println("Determine SNPs...");
			File polymorphismsFas=getPolymorphisms(alignmentFiles, core,refs.get(i));
			buildTree(polymorphismsFas);		
			if(delete){
				
				delete(alignmentFiles);
			}
		}
	}
	
	 void writeErrorFile(String errormsg){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(masterOutFolder+"/error.txt"));
			bw.write(errormsg);
			bw.close();
		}catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	private ArrayList<File> getAlignments(boolean exists){
		ArrayList<File> al=new ArrayList<File>();
		File[] files=sequenceFolder.listFiles();
		for(int i=0;i<files.length;i++){
			String name=files[i].getName();
			String suffix=aligner==soap2?".sop":PerformBowtie.samtoolsExist()?".bam":".sam";
			File alFile=null;
			String id=getId(name);
			if(hasExtension(name,fasExt)||hasExtension(name,gbkExt)){
				alFile=new File(alignmentFolder+"/"+id+"_"+arguments.get("readLength")+"fasta"+suffix);
			}else if(hasExtension(name,fastqExt)){
				//alFile=new File(alignmentFolder+"/"+id+"_"+arguments.get("readLength")+"fastq"+suffix);
				alFile=new File(alignmentFolder+"/"+id+suffix);
			}else if(hasExtension(name,gzExt)){
				alFile=new File(alignmentFolder+"/"+id+suffix);
			}
			if(alFile!=null&&!(exists ^ alFile.exists())){
				al.add(alFile);
			}
		}
		return al;
	}
	

	
	File getSequenceFile(File folder,String name){
		for(int i=0;i<fasExt.length;i++){
			File temp=new File(folder+"/"+name+"."+fasExt[i]);
			if(temp.exists())return temp;
		}
		for(int i=0;i<gbkExt.length;i++){
			File temp=new File(folder+"/"+name+"."+gbkExt[i]);
			if(temp.exists())return temp;
		}
		return null;
	}
	
	private File moveRefToCore(String ref){
		File sequenceFile=getSequenceFile(sequenceFolder,ref);
		
		File coreFolder=new File(outFolder+"/core/");
		if(!coreFolder.exists()){
			coreFolder.mkdir();
		}
		File core=new File(coreFolder+"/"+ref+".fas");
		if(!core.exists()){
			System.out.println(sequenceFile+" "+core);
			copyFasta(sequenceFile,core);
			runAlignment=true;
		}else{
			runAlignment=false;
		}
		
		return core;
	}
	
	private void copyFasta(File from,File to){
		ArrayList<Fasta> fas;
		if(hasExtension(from.getName(),gbkExt)){
			fas=new ReadGenbank(from).getSequence();
		}else{
			fas=Fasta.readFasta(from);
		}
		for(int i=0;i<fas.size();i++){
			int length=fas.get(i).getSequence().length();
			String ident=fas.get(i).getIdent().split("\\s+")[0];
			fas.get(i).setIdent(ident+" 1.."+length+" + "+ident);
		}
		Fasta.write(fas,to);
		
	}
	
	
	private ArrayList<String> getReference(){
		ArrayList<String> refIDs=new ArrayList<String>();
		File[] list=sequenceFolder.listFiles();
		for(int i=0;i<list.length;i++){
			String name=list[i].getName();
			
			if(hasExtension(name,gbkExt)){
				String id=getId(name);
				refIDs.add(id);
				return refIDs;
			}else if(hasExtension(name,fasExt)&&!(Boolean)arguments.get("genes")){
				String id=getId(name);
				refIDs.add(id);
				return refIDs;
			}
		}
		return refIDs;
	}
	
	public static boolean hasExtension(String name,String ext[]){
		for(int i=0;i<ext.length;i++){
			if(name.toLowerCase().endsWith("."+ext[i]))return true;
		}

		return false;
		
	}
	
	
	
	private ArrayList<String> getReferences(){

		ArrayList<String> refIDs=new ArrayList<String>();
		for(int i=0;i<500;i++){
			if(arguments.containsKey("ref"+i)){
				refIDs.add(arguments.get("ref"+i).toString());
			}
		}
		if(arguments.containsKey("ref")&&arguments.get("ref").toString().length()>0)refIDs.add(arguments.get("ref").toString());
		refIDs.addAll(reflist);
		return refIDs;
	}
	

	
	public void mergeAnalyses() throws RealphyException{
//		if((Boolean)arguments.get("genes")){
//			throw new RealphyException("The -genes option does currently not work in conjunction with the -merge option. I will fix this problem as soon as possible.\n");
//		}
		refs=!references?getReference():getReferences();
		mergeFolder=new File(masterOutFolder+"/merge/");
		if(!mergeFolder.exists()){
			mergeFolder.mkdir();
		}
		ArrayList<File> columnCollection=new ArrayList<File>();
		//ArrayList<Columns> totalColumns=new ArrayList<Columns>();
		for(int i=0;i<refs.size();i++){
			outFolder=new File(masterOutFolder+"/"+refs.get(i)+arguments.get("suffix"));
			clean();
			if(!outFolder.exists())outFolder.mkdir();
			String genes=(Boolean)arguments.get("genes")==false?"_NoGenes":"";
			String varOnly=(Boolean)arguments.get("varOnly")==false?"":"_noInvar";
			String base=outFolder+"/PolySeqOut"+genes+varOnly+"/";

			polymorphismsOutFolder=new File(base);
			if(!polymorphismsOutFolder.exists()){
				polymorphismsOutFolder.mkdir();
			}

			//File cutFolder=new File(sequenceFolder+"/cut/");
			alignmentFolder=!(Boolean)arguments.get("genes")?new File(outFolder+"/alignOut_NoGenes/"):new File(outFolder+"/alignOut/");
			File core=!(Boolean)arguments.get("genes")==true?moveRefToCore(refs.get(i)):determineCore(refs.get(i));
			ArrayList<File> missingAlignmentFiles=getAlignments(false);
			if(missingAlignmentFiles.size()>0){
				System.out.println("Cut fasta sequence files into fragments...");
				ArrayList<File> cuts=cutSequences();
				System.out.println("Running alignment program to align sequence fragments to reference...");
				runAlignmentProgram(core,cuts);
				//FileHandler.deleteFolder(cutFolder);
			}
			File columnsFile=new File(outFolder+"/"+refs.get(i)+"_columns.obj");
			if(!columnsFile.exists()){
				ArrayList<File> alignmentFiles=getAlignments(true);
				GetPolymorphisms gps=getColumns(alignmentFiles, core,refs.get(i));
				System.out.println("Determining polymorphic columns...");
				Clashes columns=gps.calculateColumns();
				ObjectIO.writeObject(columnsFile,columns);
				
			}
			columnCollection.add(columnsFile);
			//totalColumns.add(gps.calculateTotalColumns());
			
			
		}
		//compareColumns(columnCollection,refs);//2 specific for soap!
		//File polyCol=produceColumnStatistic(columnCollection,"polymorphicColumns.txt");
		//produceGraph(polyCol);
		//produceColumnStatistic(totalColumns,"totalColumns.txt");
		//ArrayList<Columns> cc=readColumns(columnCollection);
		File alignment=new File(mergeFolder+"/mergedAlignment.fas");
		if(columnCollection.size()>0){
			Clashes merge=new Clashes(columnCollection);
			
			merge.printAlignment(alignment);
		}
		System.out.println("Building merge tree...");
		buildTree(alignment);

	}
	
//	private void produceGraph(File polyCol){
//		String RscriptExecutable;
//		if((RscriptExecutable=getRscript())!=null){
//			File jpeg=new File(mergeFolder+"/SNPs_multi_genome.jpg");
//			File jpegTxt=new File(mergeFolder+"/SNPs_multi_genome.txt");
//			RCode rc=new RCode();
//			rc.addRCode("library(gplots)");
//			rc.addRCode("table<-read.table(\""+polyCol+"\")");
//			rc.addRCode("lt<-length(table[,1])");
//			rc.addRCode("stderr<-table[,2]");
//			rc.addRCode("miny<-min(table[,1])");
//			rc.addRCode("maxy<-max(table[,1])");
//			rc.addRCode("yl<-c(miny-1/20*miny,maxy+1/20*maxy)");
//			rc.addRCode("jpeg(\""+jpeg+"\",height=8,width=7.5,units=\"cm\",res=300)");
//			rc.addRCode("par(cex=0.5)");
//			rc.addRCode("plotCI(1:lt,table[,1],uiw=stderr,gap=0,err=\"y\",ylim=yl,ylab=\"Average number of unique polymorphic sites\", xlab=\"Number of genomes included in analysis\",axes=FALSE)");
//			rc.addRCode("axis(1,at=1:lt)");
//			rc.addRCode("axis(2)");
//			rc.addRCode("dev.off()");
//			RCaller rcall=new RCaller();
//			rcall.setRscriptExecutable(RscriptExecutable);
//			rcall.setRCode(rc);
//			//System.out.println(jpegTxt);
//			R_functions.writeCode(rc,jpegTxt);
//			rcall.runOnly();
//		}
//	}
	

	/*
	private ArrayList<Columns> readColumns(ArrayList<File> columnFiles){
		ArrayList<Columns> cc=new ArrayList<Columns>();
		for(int i=0;i<columnFiles.size();i++){
			cc.add((Columns)ObjectIO.readObject(columnFiles.get(i)));
		}
		return cc;
	}
	
	public File produceColumnStatistic(ArrayList<File> columnFiles,String name){
		File columnStats=new File(mergeFolder+"/"+name);
		try{
			ArrayList<Columns> cc=readColumns(columnFiles);
			BufferedWriter bw=new BufferedWriter(new FileWriter(columnStats));

			for(int i=1;i<cc.size()+1;i++){

				Stats s=ColumnPairStatistic.drawColumnNumbers(cc, 100, i);
				bw.write(s.getAverage()+"\t"+s.getStandardDeviation()+"\t"+s.getStandardError()+"\n");
			}
			bw.close();						
		}catch(IOException e){
			e.printStackTrace();
		}
		return columnStats;
	}
	
	public void compareColumns(ArrayList<File> cc,ArrayList<String> refs){
		File out=new File (mergeFolder+"/columnComparison.txt");
//		HashMap<String,String> genomes=readGenomes(sequenceFolder);
		for(int i=0;i<cc.size();i++){
			Columns ci=(Columns)ObjectIO.readObject(cc.get(i));
			for(int j=i+1;j<cc.size();j++){
				Columns cj=(Columns)ObjectIO.readObject(cc.get(j));
				ColumnPairStatistic cps=ColumnPairStatistic.compareColumnPair(ci,cj,refs.get(i),refs.get(j));
				cps.write(out);
//				ArrayList<String> IDs=cps.getDistinct();
//				if((Boolean)arguments.get("evaluateOverlaps")){
//					EvaluateRegions er=new EvaluateRegions(IDs,genomes,maxPol,sequenceFolder,outFolder,flank+((Integer)arguments.get("perBaseCov")/2+2));
//					er.write(out);
//				}
			}
		}
	}*/
	
//	private static HashMap<String,String> readGenomes(File sequenceFolder){
//		HashMap<String,String> genomes=new HashMap<String, String>();
//		File[] list=sequenceFolder.listFiles();
//		for(int i=0;i<list.length;i++){
//			String name=list[i].getName();
//			if(name.endsWith(".fas")){
//				genomes.putAll(Fasta.fasToHash(Fasta.readFasta(list[i])));
//			}
//		}
//		return genomes;
//	}

//	//determine core genome
//	public File determineCore(){
//		File outFolderCore=new File(outFolder+"/core/");
//		if(!outFolderCore.exists()){
//			if(outFolderCore.mkdir()==false){
//				System.err.println(outFolderCore+" could not be created.");
//				System.exit(-1);
//			}
//			runSoap=true;
//		}
//		File core=getCoreFile(outFolderCore);
//		if(core==null){
//			String[] dcgArgs={sequenceFolder.toString(),outFolderCore.toString(),arguments.get("perBaseScore").toString(),arguments.get("readLength").toString(),"blastn",arguments.get("coreGenomes").toString(),1+"",path.get("BLASTBINPATH"),arguments.get("ref").toString()};
//			core=DetermineCoreGenome.runDetermineCoreGenome(dcgArgs).get(0);
//			runSoap=true;
//		}
//		return core;
//	}
	public File determineCore(String ref){
		File outFolderCore=new File(outFolder+"/core/");
		if(!outFolderCore.exists()){
			outFolderCore.mkdir();
			runAlignment=true;
		}
		File core=clean?null:getCoreFile(outFolderCore);
		if(core==null){
			String[] dcgArgs={sequenceFolder.toString(),outFolderCore.toString(),arguments.get("perBaseScore").toString(),arguments.get("readLength").toString(),"blastn",arguments.get("coreGenomes").toString(),1+"",path.get("BLAST"),path.get("BLASTBUILDER"),ref};
			core=DetermineCoreGenome.runDetermineCoreGenome(dcgArgs).get(0);
			runAlignment=true;
		}
		return core;
	}
	public ArrayList<File> cutSequences(){
		//cut up sequences
		String[] cuArgs={sequenceFolder.toString(),arguments.get("readLength").toString(),cutFolder.toString(),new Boolean(clean).toString()};
		ArrayList<File> cutSequences=CutUpSequences.runCutUpSequences(cuArgs);
		
		return cutSequences;
	}
	
	public ArrayList<File> runAlignmentProgram(File core,ArrayList<File> cutSequences){
		//soap sequences

		if(!alignmentFolder.exists()){
			alignmentFolder.mkdir();
		}
		try{
			ArrayList<File> alignmentFiles=aligner==soap2?PerformSoap.runMultipleSoaps(core,cutSequences,alignmentFolder,path.get("SOAP"),path.get("SOAPBUILDER"),runAlignment):PerformBowtie.runMultiple(core,cutSequences,alignmentFolder,path.get("BOWTIE2"),path.get("BOWTIE2BUILDER"),runAlignment,seedLength,bowtieOptions);
			return alignmentFiles;
		}catch(RealphyException e){
			e.printStackTrace();
			writeErrorFile(e.getMessage());
			System.exit(-1);
			return null;

		}
		
	}

	public void delete(ArrayList<File> files){
		File folder=files.get(0).getParentFile();
		for(int i=0;i<files.size();i++){
			if(!files.get(i).delete())System.err.println("File "+files.get(i)+" could not be deleted.");
		}
		folder.delete();
	}
	public GetPolymorphisms getColumns(ArrayList<File> alignmentFiles,File core,String ref){
		//determine polymorphisms
		System.out.println("Determine polymorphisms...");
		try{
			GetPolymorphisms GPS=new GetPolymorphismWithGaps(alignmentFiles,refs, core, flank, (Integer)arguments.get("quality"), (Double)arguments.get("polyThreshold"), (Double)arguments.get("fractionCov"), perBaseCov,(Boolean)arguments.get("merge"),!(Boolean)arguments.get("genes"),gapThreshold,!(Boolean)arguments.get("varOnly"),polymorphismsOutFolder,ref);

			File polymorphismsFas=GPS.writeSequences(); 
			System.out.println("Building tree...");
			buildTree(polymorphismsFas);
			return GPS;
		}catch(RealphyException e){
			e.printStackTrace();
			writeErrorFile(e.getMessage());
			System.exit(-1);
			return null;
		}
	} 
	
	public File getPolymorphisms(ArrayList<File> alignmentFiles,File core,String ref){
		//determine polymorphisms
		try{	
			GetPolymorphisms GPS=new GetPolymorphismWithGaps(alignmentFiles, refs,core, flank, (Integer)arguments.get("quality"), (Double)arguments.get("polyThreshold"), (Double)arguments.get("fractionCov"), perBaseCov,(Boolean)arguments.get("merge"),!(Boolean)arguments.get("genes"),gapThreshold,!(Boolean)arguments.get("varOnly"),polymorphismsOutFolder,ref);

			File polymorphismsFas=GPS.writeSequences();
			return polymorphismsFas;
		}catch(RealphyException e){
			e.printStackTrace();
			writeErrorFile(e.getMessage());
			System.exit(-1);
			return null;
		}
	}
	
	public void buildTree(File polymorphismsFas){
		
		//build tree
			//to phylip format
		String root=((String)arguments.get("root")).split(",")[0];	
		File polymorphismsMove=new File(polymorphismsFas.getParent()+"/polymorphisms_move.fas");
		int seqLength=MoveSeqToTop.move(root,polymorphismsFas,polymorphismsMove);
		polymorphismsFas.delete();
		File phylip=new File(polymorphismsMove.getParent()+"/polymorphisms_move.phy");
		int treeBuilder=(Integer)arguments.get("treeBuilder");
		ArrayList<Fasta> fas=Fasta.readFasta(polymorphismsMove);
		if(treeBuilder>0 && fas.size()>3){
			System.out.println("Build tree...");
			File tree=null;
			if(treeBuilder==1){
				Fasta.writePhylip(fas, phylip, 10);
				RunTreePrograms.runTreePuzzle(phylip,path.get("TREEPUZZLE"));
				tree=new File(phylip+".tree");
			}
			else if(treeBuilder==2){
				Fasta.writePhylip(fas, phylip, 100);
				int rand=(int)(Math.random()*1000000);
				RunTreePrograms.runRAxML(phylip.getAbsoluteFile(), new File(path.get("RAXML")),seqLength,"raxml",(String)arguments.get("root"),!(Boolean)arguments.get("genes"),treeOptions,rand);
				tree=new File(phylip.getParent()+"/RAxML_bestTree.raxml");
			}else if(treeBuilder==3){
				Fasta.writePhylip(fas, phylip, 10);
				tree=RunTreePrograms.runMaxPars(phylip.getAbsoluteFile(), new File(path.get("MaxPars")));
			}else if(treeBuilder==4){
				Fasta.writePhylip(fas, phylip, 100);
				int rand=(int)(Math.random()*1000000);
				RunTreePrograms.runPhyML(phylip.getAbsoluteFile(), new File(path.get("PhyML")),treeOptions,rand);
				tree=new File(phylip+"_phyml_tree.txt");
			}
			if(tree!=null)printTree(tree);
		}else if(treeBuilder>0){
			String errormsg="TOO FEW SPECIES! For tree reconstruction at least four different isolates are required.";
			System.err.println("TOO FEW SPECIES! For tree reconstruction at least four different isolates are required.");
            writeErrorFile(errormsg);
		}
	}
		
	private void printTree(File tree){
		String RscriptExecutable=path.get("Rscript");
		if(RscriptExecutable==null){
			System.err.println("Cannot find Rscript executable. Please set \"Rscript\" in config.txt.");
			return;
		}
		File jpeg=new File(tree.getParent()+"/tree.jpg");
		RCode rc=R_functions.plot_InitJpg(jpeg);
		R_functions.plotTree(rc,tree);
		R_functions.runRCode(rc,new File(RscriptExecutable));
		//R_functions.writeCode(rc,new File(jpeg+".txt"));

	}
	
	private File getCoreFile(File folder){
		File[] list=folder.listFiles();
			for(int i=0;i<list.length;i++){
				String file=list[i].toString();
				if(file.endsWith("Flank"+flank+".fas")){
					return list[i];
				}
			}
		
		return null;
		
	}
	

	
//	private String adjustRootLengths(String root){
//		String[] split=root.split(",");
//		StringBuffer newRoot=new StringBuffer(cutString(split[0]));
//		for(int i=1;i<split.length;i++){
//			newRoot.append(","+cutString(split[i]));
//		}
//		return newRoot.toString();
//	}
//	
//	private String cutString(String cut){
//		if(cut.length()>10){
//			return cut.substring(0,10);
//		}else return cut;
//	}
	

	
	
	
	
	private void initDefault(){
		path.put("SOAP", "soap");
		path.put("TREEPUZZLE", "puzzle");
		path.put("RAXML", "raxmlHPC-SSE3");
		path.put("JAVA", "java");
		path.put("BLAST","blastall");
		path.put("BLASTBUILDER","formatdb");
		path.put("SOAPBUILDER","2bwt-builder");
		path.put("BOWTIE2","bowtie2");
		path.put("BOWTIE2BUILDER","bowtie2-build");
		path.put("MaxPars","dnapars");
		path.put("PhyML","phyml");
		arguments.put("readLength", new Integer(50));
		arguments.put("coreGenomes", new Integer(1));
		arguments.put("perBaseScore", new Double(1));
		arguments.put("quality",new Integer(20));
		arguments.put("polyThreshold",new Double(0.95));
		arguments.put("fractionCov",new Double(0));
		arguments.put("perBaseCov",new Integer(10));
		arguments.put("root","");
		arguments.put("evaluateOverlaps",new Boolean(false));
		arguments.put("merge",new Boolean(false));
		arguments.put("genes", new Boolean(false));
		arguments.put("gapThreshold", new Double(0));
		arguments.put("treeBuilder", new Integer(4));
		arguments.put("clean", new Boolean(false));
		arguments.put("quiet", new Boolean(false));
		arguments.put("varOnly", new Boolean(false));
		//arguments.put("covWindow", 0);
		arguments.put("aligner", bowtie2);
		arguments.put("seedLength", 22);
		arguments.put("suffix","");
		arguments.put("gaps", false);


	}
	

	
	private void setVariables(String[] args){
		sequenceFolder=new File(args[0]);
		masterOutFolder=new File(args[1]);
		cutFolder=new File(sequenceFolder+"/cut/");
		int i=2;
		while(i<args.length){
			int add=0;	
			if(!args[i].startsWith("-")){
				System.err.println("Arguments have to start with \"-\". Violating argument: "+args[i]+".");
				System.exit(-1);
			}else{
				String arg=args[i].substring(1);
				if(novalueArgs.containsKey(arg)){
					add=1;
					Object o=new Object();
					if((o=checkValue(arg,""))!=null){
						arguments.put(arg, o);

					}else{
						System.err.println("Invalid value \""+args[i+1]+"\" for argument -"+arg+".");
						System.exit(-1);
					}
				}else{
					add=2;
					Object o=new Object();
					if((o=checkValue(arg,args[i+1]))!=null){
						arguments.put(arg, o);

					}else{
						System.err.println("Invalid value \""+args[i+1]+"\" for argument -"+arg+".");
						System.exit(-1);
					}
				}
			}
			i+=add;
		}
		checkReference();
		config=config==null?new File(masterOutFolder+"/config.txt"):config;
		if(config.exists()){
			readConfig(config);
		}else{
			System.err.println("Cannot find "+config+"!");
			System.exit(-1);
		}

	}
	
	private void checkReference(){
		for(int i=0;i<502;i++){
			String ref="";
			if(arguments.containsKey("ref"+i)||(i==501&&arguments.containsKey("ref"))){
				ref=i==501?(String)arguments.get("ref"):(String)arguments.get("ref"+i);
				File seqFile=getSequenceFile(sequenceFolder, ref);
				if(seqFile!=null){
					boolean fas=hasExtension(seqFile.getName(), fasExt);
					boolean gbk=hasExtension(seqFile.getName(), gbkExt);
					if(!gbk&&(Boolean)arguments.get("genes")){
						System.err.print("File called ");
						for(int j=0;j<gbkExt.length-1;j++){
							System.err.println(ref+"."+gbkExt[j]+" or ");
						}
						System.err.println(ref+"."+gbkExt[gbkExt.length-1]+" is required as reference with CDS information if genes is set.");
						System.exit(-1);
					}else
						if(!fas&&!gbk && !(Boolean)arguments.get("genes")){
							System.err.print("File called ");
							for(int j=0;j<fasExt.length-1;j++){
								System.err.print(ref+"."+fasExt[j]+" or ");
							}
							System.err.println(ref+"."+fasExt[fasExt.length-1]+" is required as reference.");
							System.exit(-1);
						}
				}else{
					System.err.println("Cannot find reference file with the prefix: "+ref+".");
				}
			}
		}
	}
	
	public void readConfig(File in){
		try{
			BufferedReader br=new BufferedReader(new FileReader(in));
			String line="";
			while((line=br.readLine())!=null){
				String[] split=line.split("\\s+");
				if(line.length()<4)continue;
				path.put(split[0],split[1]);
				File temp=new File(split[1]);
				if(!temp.exists()){
					System.err.println("Cannot find file "+temp+"!");
					if(aligner==soap2){
						if(split[0].equals("SOAP")||split[0].equals("SOAPBUILDER")){
							System.err.println("Cannot find "+split[1]+". The executable for "+split[0]+" is essential! Exiting.");
							System.exit(-1);
						}
					}else if(aligner==bowtie2){
						if(split[0].equals("BOWTIE2")||split[0].equals("BOWTIE2BUILDER")){
							System.err.println("Cannot find "+split[1]+". The executable for "+split[0]+" is essential! Exiting.");
							System.exit(-1);
						}
					}
				}
			}
			
			br.close();
		}catch(IOException e){
			e.printStackTrace();
		}
	}
	
	public Object checkValue(String arg,String value){
		if(arg.equals("readLength")){
			return checkInteger(arg,value,30,Integer.MAX_VALUE);
		}else if(arg.equals("coreGenomes")){
			return checkInteger(arg,value,1,Integer.MAX_VALUE);
		}else if(arg.equals("perBaseScore")){
			return checkDouble(arg,value,0,2);
		}else if(arg.equals("quality")){
			return checkInteger(arg,value,0,41);
		}else if(arg.equals("polyThreshold")){
			return checkDouble(arg,value,0,1);
		}else if(arg.equals("fractionCov")){
			return checkDouble(arg,value,0,1);
		}else if(arg.equals("gapThreshold")){
			gapThreshold=checkDouble(arg,value,0,1.0);
			return gapThreshold;
		}else if(arg.equals("gaps")||arg.equals("g")){
			gaps=true;
			if(gaps&&gapThreshold==0){
				gapThreshold=1.1;
			}
			
			return true;

		}else if(arg.equals("perBaseCov")){
			Integer i;
			if((i=checkInteger(arg,value,1,Integer.MAX_VALUE))!=null){
				perBaseCov=i;
				return i;
			}else{
				return null;
			}
		}else if(arg.equals("treeBuilder")){
			return checkInteger(arg,value,0,4);
		}else if(arg.equals("aligner")){
			return checkInteger(arg,value,soap2,bowtie2);
		}else if(arg.equals("covWindow")){
			return checkInteger(arg,value,0,Integer.MAX_VALUE);
		}else if(arg.equals("seedLength")){ 
			return checkInteger(arg,value,4,32);
		}else if(arg.equals("ref")){
			references=true;
			return checkRef(value);

		}else if(arg.equals("ref-list")){
			references=true;
			return checkReflist(value);

		}else if(arg.equals("root")){
			return checkRoot(value);

		}else if(arg.startsWith("ref")){
			references=true;
			return checkRef(value);

		}else if(arg.startsWith("config")){
			return checkConfig(value);

		}else if(arg.equals("merge")||arg.equals("m")){
			arguments.put("merge", true);

			return true;

		}else if(arg.equals("genes")){
			arguments.put("genes", true);

			return true;

		}else if(arg.startsWith("clean")||arg.equals("c")){
			arguments.put("clean", true);

			return true;

		}else if(arg.equals("quiet")||arg.equals("q")){
			arguments.put("quiet", true);

			return true;

		}else if(arg.equals("varOnly")||arg.equals("v")){
			arguments.put("varOnly", true);

			return true;

		}else if(arg.equals("suffix")){
			return value;

		}else if(arg.equals("d")||arg.equals("delete")){
			arguments.put("d", true);
			return true;

		}else if(arg.equals("version")){
			printVersion();
			return null;
		}else if(arg.equals("h")||arg.equals("help")){
			printHelp();
			System.exit(0);
			return null;
		}else if(arg.equals("treeOptions")){
			treeOptions=checkFile(value);
			return treeOptions;
		}else if(arg.equals("bowtieOptions")){
			bowtieOptions=checkFile(value);
			return bowtieOptions;
		}else{
			System.err.println("Do not recognize "+arg+".");
			return null;
		}
	}
	
	public Double checkDouble(String arg,String value,double l,double u){
		try{
			double x=Double.parseDouble(value);
			if(x<l){
				System.err.println(arg+" needs to be greater than "+l+".");
				return null;
			}else if(x>u){
				System.err.println(arg+" needs to be smaller than "+u+".");
				return null;
			}else{
				return x;
			}
		}catch(NumberFormatException e){
			System.err.println(value+" is not a valid double.");
			return null;
		}
		
	}
	
	public Boolean checkBoolean(String bool){
			return Boolean.parseBoolean(bool);
		
		
	}
	
	public String checkRoot(String root){
		String temp[]=root.split(",");
		for(int i=0;i<temp.length;i++){
			
			File tempFile=getSequenceFile(sequenceFolder,temp[i]);
			if(tempFile==null){
				System.err.println("No file with prefix "+temp[i]+" found.");
				return null;
			}
		}
		return root;
	}
	
	public String checkRef(String ref){
		File temp=getSequenceFile(sequenceFolder, ref);
		if(temp!=null){
			return ref;
		}else{
			System.err.println("No file with prefix "+ref+" exists.");
			return null;
		}
	}
	public String checkConfig(String configString){
		config=new File(configString);
		boolean conf=config.exists();
		if(conf){
			return configString;
		}else{
			System.err.println("No file called "+configString+" exists.");
			return null;
		}
	}
	public Boolean checkReflist(String refListFile){
		reflist=readReferenceList(new File(refListFile));
		for(int i=0;i<reflist.size();i++){
			checkRef(reflist.get(i));
		}
		return true;
	}
	
	private ArrayList<String> readReferenceList(File in){
		ArrayList<String> refs=new ArrayList<String>();
		try{
			BufferedReader br=new BufferedReader(new FileReader(in));
			String line;
			while((line=br.readLine())!=null){
				refs.add(getId(line));
			}
			
			
			br.close();
		}catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
		return refs;
	}
	

	
	public File checkFile(String file){
		File in=new File(file);
		if(in.exists()){
			return in;
		}else{
			throw new RuntimeException("File "+file+" does not exist.");
		}
	}
	
	public Integer checkInteger(String arg,String value,int l,int u){
		try{
			int x=Integer.parseInt(value);
			if(x<l){
				System.err.println(arg+" needs to be greater/equal than "+l+".");
				return null;
			}else if(x>u){
				System.err.println(arg+" needs to be smaller/equal than "+u+".");
				return null;
			}else{
				return x;
			}
		}catch(NumberFormatException e){
			System.err.println(value+" is not a valid integer.");
			return null;
		}
		
	}
	
}
