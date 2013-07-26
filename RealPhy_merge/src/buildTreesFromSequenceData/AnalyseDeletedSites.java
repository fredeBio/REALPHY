package buildTreesFromSequenceData;

import java.io.*;
import java.util.*;

import rcaller.RCode;

import util.R_functions;

public class AnalyseDeletedSites {
	public static void main(String args[]){
		File folder=new File(args[0]);
		int covThreshold=Integer.parseInt(args[1]);
		File scriptPath=new File(args[2]);
		double polT=Double.parseDouble(args[3]);
		AnalyseDeletedSites ads=new AnalyseDeletedSites(covThreshold, polT, folder);
		
		ads.plotHistograms(new File(folder+"/deletedHistograms.jpg"),scriptPath);
		ads.writeDeletedAbovePolT(new File(folder+"/deletedTop"+polT+".txt"));
		
	}
	
	public AnalyseDeletedSites(int covThreshold,double polT,File folder){
		this.covThreshold=covThreshold;
		this.polT=polT;
		readDeletedSites(folder);
	}
	
	//minimum coverage for sites to be included in analysis
	int covThreshold=10;
	
	//only include sites that have at least a proportion of polyT polymophisms of coverage in the second part of the analysis
	double polT=0.9;
	
	ArrayList<String> names=new ArrayList<String>();
	
	//list of different querygenomes<list of (major polymorphisms/coverage)>
	ArrayList<ArrayList<Double>> data=new ArrayList<ArrayList<Double>>();

	public void readDeletedSites(File folder){
		File[] files=folder.listFiles();
		for(int i=0;i<files.length;i++){
			String name=files[i].getName();
			if(name.endsWith("deletedSites.txt")){
				ArrayList<Double> ratios=readRatio(files[i],covThreshold);
				if(ratios.size()>0){
					data.add(ratios);
					names.add(name);
				}
			}
		}
	}
	
	public void writeDeletedAbovePolT(File out){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			for(int i=0;i<names.size();i++){
				int count=0;
				for(int j=0;j<data.get(i).size();j++){
					if(data.get(i).get(j)>=polT){
						count++;
					}
				}
				bw.write(names.get(i)+"\t"+count+"\n");
			}
			bw.close();
		}catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	
	public  void plotHistograms(File out,File scriptPath){
		
		int col=6;
		int row=data.size()/col+1;
		RCode rc=R_functions.plot_InitJpg(out, row, col);
		
		for(int i=0;i<data.size();i++){
				R_functions.makeHistogram(rc,data.get(i),0,1,10,names.get(i));
		}
		R_functions.runRCode(rc, scriptPath);
		R_functions.writeRCode(rc, new File(out+".txt"));
	}
	
	public static ArrayList<Double> readRatio(File del,int covThreshold){
		ArrayList<Double> ratios=new ArrayList<Double>();
		try{
			
			BufferedReader br=new BufferedReader(new FileReader(del));
			String line;
			//System.out.println(del);

			while((line=br.readLine())!=null){
				String split[]=line.split("\\s+");
				//if(split.length>4){
					double cov=Double.parseDouble(split[3]);
					if(cov>=covThreshold&& !split[5].equals("NaN")){
						double ratio=Double.parseDouble(split[5]);
						ratios.add(ratio);
					}
				//}
				
			}
			br.close();
		}catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
		return ratios;
	}
}
