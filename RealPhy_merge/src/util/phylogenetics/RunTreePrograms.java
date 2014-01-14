package util.phylogenetics;

import java.io.*;
import java.util.*;

public class RunTreePrograms {
	public static void main(String args[]){
		runProgram("/home/frederic/Programs/RAxML/raxmlHPC-SSE3 -s /home/frederic/test/PolySeqOut_NoGenes_noInvar/polymorphisms_move.phy -w /home/frederic/test/PolySeqOut_NoGenes_noInvar -m GTRGAMMA -p 1234 -n inc1  -o S11 -f e -t /home/frederic/test/PolySeqOut_NoGenes_noInvar/inc1.tree", "", new File("/home/frederic/test/"), 5000);
	}
	public static File runMaxPars(File alignmentPhy,File maxParsPath){
		//System.out.println();
		File out=new File(alignmentPhy.toString().split("\\.")[0]+".tree");
		File outFolder=alignmentPhy.getParentFile();
		String command=maxParsPath.toString();
		RunTreePrograms.deleteConsenseFiles(outFolder);
		String input=alignmentPhy+"\nY\n";
		//System.out.println(command+" "+input+" "+outFolder);
		RunTreePrograms.runProgram(command, input, outFolder);
		File outtree=new File(outFolder+"/outtree");
		outtree.renameTo(out);
		return out;
	}

	
	public static void runPhyML(File in,File phymlPath,File parameters,int seed){
		File outFolder=new File(in.getParent());
		if(!phymlPath.exists()){
			System.err.println("Cannot find phyml executable. Not building tree. Exiting.");
			System.exit(-1);
		}		

		HashMap<String,Boolean> paraHM=getParameters(parameters);
		String paraLine=getParametersLine(parameters);
		String model=paraHM.containsKey("-m")?"":" -m GTR ";
		String seedString=paraHM.containsKey("--r_seed")?"":" --r_seed "+seed+" ";
		String bootstrap=paraHM.containsKey("-b")?"":" -b 0 ";
		String PhyMLcom=phymlPath+" -i "+in+model+seedString+" "+paraLine+bootstrap;
		
		runProgram(PhyMLcom,"",outFolder);
}
	
	public static void runRAxML(File in,File RAxMLPath,int seqLength,String suffix,String root,boolean noGenes,File parameters,int seed){
			File outFolder=new File(in.getParent());
			if(!RAxMLPath.exists()){
				System.err.println("Cannot find raxml executable. Not building tree. Exiting.");
				System.exit(-1);
			}
		
			String modelFile="";			
			deleteAllRAxMLFiles(outFolder,suffix);
			
			if(root.length()>0){
				root=" -o "+root;
			}

			if(noGenes==false){
				File model=generateModelFile(seqLength,outFolder);
				modelFile=" -q "+model.toString();
			}
			HashMap<String,Boolean> paraHM=getParameters(parameters);
			String paraLine=getParametersLine(parameters);
			String model=paraHM.containsKey("-m")?"":" -m GTRGAMMA ";
			String seedString=paraHM.containsKey("-p")?"":" -p "+seed+" ";
			root=paraHM.containsKey("-o")?"":root;
			String raxMLcom=RAxMLPath+" -s "+in+" -w "+outFolder+model+seedString+" -n "+suffix+" "+root+" "+modelFile+" "+paraLine;
			
			runProgram(raxMLcom,"",outFolder);
	}

	public static HashMap<String,Boolean> getParameters(File in){
		HashMap<String,Boolean> para=new HashMap<String, Boolean>();
		try{
			if(in!=null&&in.exists()){
				BufferedReader br=new BufferedReader(new FileReader(in));
				String line=br.readLine();
				String split[]=line.split("\\s+");
				for(int i=0;i<split.length;i++){
					para.put(split[i], true);
				}
				br.close();
			}
		}catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
		return para;
	}
	
	public static String getParametersLine(File in){
		String para="";
		try{
			if(in!=null&&in.exists()){
				BufferedReader br=new BufferedReader(new FileReader(in));
				String line=br.readLine();
				para=line;
				br.close();
			}
		}catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
		return para;
	}
	
	public static File runRAxMLGivenTopo(File topo,File in,File RAxMLPath,String suffix,String root,int seed,boolean delete,long time){
		File outFolder=new File(in.getParent());
		if(!RAxMLPath.exists()){
			System.err.println("Cannot find raxml executable. Not building tree. Exiting.");
			System.exit(-1);
		}
		File out=new File(outFolder+"/RAxML_result."+suffix);
		if(delete||!out.exists()){
			deleteAllRAxMLFiles(outFolder,suffix);

			if(root.length()>0){
				root=" -o "+root;
			}
			String raxMLcom=RAxMLPath+" -s "+in+" -w "+outFolder+" -m GTRGAMMA -p "+seed+" -n "+suffix+" "+root+" -f e -t "+topo;
			//System.out.println(raxMLcom);

			if(runProgram(raxMLcom,"",outFolder,time)!=0){
				return null;
			}
		}
		return out;
	}

	public static File runRAxMLSiteLikelihood(File tree,File in,File RAxMLPath,String suffix,String root,int seed,boolean delete){
		File outFolder=new File(in.getParent());
		if(!RAxMLPath.exists()){
			System.err.println("Cannot find raxml executable. Not building tree. Exiting.");
			System.exit(-1);
		}
		File out=new File(outFolder+"/RAxML_perSiteLLs."+suffix);
		if(delete||!out.exists()){
			deleteAllRAxMLFiles(outFolder,suffix);
			if(root.length()>0){
				root=" -o "+root;
			}

			String raxMLcom=RAxMLPath+" -s "+in+" -w "+outFolder+" -m GTRGAMMA -p "+seed+" -n "+suffix+" -f g -z "+tree;

			runProgram(raxMLcom,"",outFolder);
		}
		return out;
	}
	public static int runProgram(String command,String input,File dir){
		return runProgram(command,input,dir,null,null,-1);


	}
	public static int runProgram(String command,File input,File dir){
		return runProgram(command,"",dir,null,input,-1);


	}
	public static int runProgram(String command,String input,File dir,long time){
		return runProgram(command,input,dir,null,null,time);


	}
	public static int runProgram(String command,String input,File dir,File standardOut){
		return runProgram(command,input,dir,standardOut,null,-1);


	}
	public static int runProgram(String command,String input,File dir,File standardOut,File standardIn,long time){
		Timer timer = null;
		Process p=null;
		Boolean[] timeOut=new Boolean[]{false};

		try{

			Runtime runtime=Runtime.getRuntime();
			p=runtime.exec(command,new String[0],dir);

			OutputStream stdin=p.getOutputStream();
			if(standardIn==null){
				stdin.write(input.getBytes());
				stdin.flush();
				stdin.close();
			}else{
				BufferedReader br=new BufferedReader(new FileReader(standardIn));
				BufferedWriter bw=new BufferedWriter(new OutputStreamWriter(stdin));
				String line;
				while((line=br.readLine())!=null){
					bw.write(line+"\n");
					//System.out.println(line);
				}
				br.close();
				bw.close();
				
			}
			BufferedReader bri = new BufferedReader(new InputStreamReader(p.getInputStream()));
			
			String line="";
            if(time>-1){
            	timer = new Timer(true);
            	InterruptTimerTask interrupter = new InterruptTimerTask(p,timeOut);
            	timer.schedule(interrupter, time);
            }
			BufferedReader br=new BufferedReader(new InputStreamReader(p.getErrorStream()));

			while((line=br.readLine())!=null){
				System.err.println(line);
			}
            if(standardOut==null){
            	System.out.println(command);
            	while ((line = bri.readLine()) != null) {
            		System.out.println(line);
            	}
            }else{
            	BufferedWriter bw=new BufferedWriter(new FileWriter(standardOut));
            	while ((line = bri.readLine()) != null) {
            		bw.write(line+"\n");
            	}
            	bw.close();
            }
			bri.close();
			if(p.waitFor()!=0){
				System.err.println(command);
				

				while((line=br.readLine())!=null){
					System.err.println(line);
				}
				br.close();
				
				System.err.println(command+" was not successful!");
				System.exit(-1);
			}
		}catch(InterruptedException e){
			p.destroy();
			e.printStackTrace();
			return -1;
		}catch(IOException e){
			if(!timeOut[0]){
				e.printStackTrace();
				System.exit(-1);
			}else{
				System.err.println("Program: "+command+" timed out!!!");
			}
			return -1;
		}
        finally
        {
            if(timer!=null)timer.cancel();     // If the process returns within the timeout period, we have to stop the interrupter
                                // so that it does not unexpectedly interrupt some other code later.

            Thread.interrupted();   // We need to clear the interrupt flag on the current thread just in case
                                    // interrupter executed after waitFor had already returned but before timer.cancel
                                    // took effect.
                                    //
                                    // Oh, and there's also Sun bug 6420270 to worry about here.
        }
		return 0;
	}
	public static File runConsense(File in,File consensePath){

		if(!consensePath.exists()){
			System.err.println("Cannot find "+consensePath+" executable. Not building tree. Exiting.");
			System.exit(-1);
		}


		String consensecom=consensePath.toString();
		deleteConsenseFiles(in.getParentFile());
		if(!in.getName().equals("intree"))runProgram(consensecom,in.toString()+"\ny\r\n",in.getParentFile());
		else{
			System.err.println("Please do not choose \"intree\" as file name!");
			
		}
		return new File(in.getParentFile()+"/outfile");
	}
	
	
	public static void deleteConsenseFiles(File outFolder){
		File outfile=new File(outFolder+"/outfile");
		File outtree=new File(outFolder+"/outtree");
		File intree=new File(outFolder+"/intree");
		if(outfile.exists())outfile.delete();
		if(outtree.exists())outtree.delete();
		if(intree.exists())intree.delete();
	}
	
	public static void runPhyml(File in,String PhymlPath){
		
		File phymlExec=new File(PhymlPath);
		if(!phymlExec.exists()){
			System.err.println("Cannot find Phyml executable. Not building tree. Exiting.");
			System.exit(-1);
		}
		String phymlcom=PhymlPath+" -i "+in+" -m F81 -c 1 -f e";
		runProgram(phymlcom,"",in.getParentFile());
	}
	
	public static void runDNAML(File in,String dnamlPath){

		File phymlExec=new File(dnamlPath);
		if(!phymlExec.exists()){
			System.err.println("Cannot find Phyml executable. Not building tree. Exiting.");
			System.exit(-1);
		}
		String phymlcom=dnamlPath;
		runProgram(phymlcom,in.toString()+"",in.getParentFile());
	}

	public static void runTreePuzzle(File in,String treePuzzlePath){
			
			File puzzleExec=new File(treePuzzlePath);
			if(!puzzleExec.exists()){
				System.err.println("Cannot find TREE-PUZZLE executable. Not building tree. Exiting.");
				System.exit(-1);
			}
			String treepuzzlecom=treePuzzlePath+" "+in;
			runProgram(treepuzzlecom,"w\r\ny\r\n",in.getParentFile());
	}
	private static void deleteAllRAxMLFiles(File folder,String suffix){
		File[] list=folder.listFiles();
		for(int i=0;i<list.length;i++){
			if(list[i].getName().endsWith(suffix)){
				list[i].delete();
			}
		}
	}
	
	private static File generateModelFile(int seqLength,File outFolder){
		try{
			File modelOut=new File(outFolder+"/model.txt");
			BufferedWriter bw=new BufferedWriter(new FileWriter(modelOut));
			bw.write("DNA, codon1 = 1-"+seqLength+"\\3\n"+
					"DNA, codon2 = 2-"+seqLength+"\\3\n"+
					"DNA, codon3 = 3-"+seqLength+"\\3\n");
			
			bw.close();
			return modelOut;
		}catch(IOException e){
			e.printStackTrace();
			return null;
		}
	}
}
/**
 * Just a simple TimerTask that interrupts the specified thread when run.
 */
class InterruptTimerTask
        extends TimerTask
{

    private Process p;
    private Boolean[] TimeOut;
    public InterruptTimerTask(Process p,Boolean[] timeOut)
    {
    	TimeOut=timeOut;
        this.p = p;
    }

    public void run()
    {		TimeOut[0]=true;
    		p.destroy();
    }

}