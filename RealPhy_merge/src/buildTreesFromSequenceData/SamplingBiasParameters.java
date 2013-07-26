package buildTreesFromSequenceData;

import java.io.*;
import java.util.ArrayList;


public final class SamplingBiasParameters {
	ArrayList<Integer> seqLength=new ArrayList<Integer>();
	ArrayList<Double> A=new ArrayList<Double>();
	ArrayList<Double> B=new ArrayList<Double>();
	public SamplingBiasParameters(){
		readParameters();
	}
	
	public void readParameters(){
		try{
			BufferedReader br=new BufferedReader(new InputStreamReader(this.getClass().getResourceAsStream("/resources/parameters.txt")));
			String line;
			while((line=br.readLine())!=null){
				String split[]=line.split("\\s+");
				seqLength.add(Integer.parseInt(split[0]));
				A.add(Double.parseDouble(split[1]));
				B.add(Double.parseDouble(split[2]));
			}
			br.close();
		}catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	public  double getA(int seqlength){
		for(int i=0;i<seqLength.size();i++){
			if(seqlength<=seqLength.get(i)){
				return A.get(i);
			}
		}
		return A.get(A.size()-1);
	}
	
	public double getB(int seqlength){
		for(int i=0;i<seqLength.size();i++){
			if(seqlength<=seqLength.get(i)){
				return B.get(i);
			}
		}
		return B.get(B.size()-1);
	}

}
