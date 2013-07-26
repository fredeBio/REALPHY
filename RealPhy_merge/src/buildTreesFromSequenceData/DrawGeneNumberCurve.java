package buildTreesFromSequenceData;

import java.io.*;
import java.util.*;

import statistics.Stats;

public class DrawGeneNumberCurve {
	public static void main(String args[]){
		File groupsFile=new File(args[0]);
		int draws=Integer.parseInt(args[1]);
		String id=groupsFile.getName().split("\\.")[0];
		File absoluteGeneNumbers=new File(groupsFile.getParent()+"/"+id+"_geneNumbers_"+draws+".plot");
		File numberOfGroups=new File(groupsFile.getParent()+"/"+id+"_groupNumbers_"+draws+".plot");
		ArrayList<String> strains=new ArrayList<String>();
		HashMap<String,Integer> geneNumbers=new HashMap<String,Integer>();
		HashMap<String,HashMap<Integer,Boolean>> groupNumbers=new HashMap<String,HashMap<Integer,Boolean>>();

		getGeneNumbers(groupsFile,strains,geneNumbers);
		strains=new ArrayList<String>();
		getGroups(groupsFile,strains,groupNumbers);
		System.out.println("Group file analysed.");
		calculateGeneNumbers(strains,geneNumbers,absoluteGeneNumbers,draws);
		calculateGroupNumbers(strains,groupNumbers,numberOfGroups,draws);

		
	}
	
	public static void calculateGroupNumbers(ArrayList<String> strains,HashMap<String,HashMap<Integer,Boolean>> groupNumbers,File out,int draws){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			for(int i=0;i<strains.size();i++){
				Stats s=drawGroupNumbers(strains,groupNumbers,draws,i+1);
				bw.write(s.getAverage()+"\t"+s.getStandardDeviation()+"\n");
				System.out.println("Average of "+(i+1)+" genomes: "+s.getAverage()+" +- "+s.getStandardDeviation());
			}
			bw.close();
		}catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	public static void calculateGeneNumbers(ArrayList<String> strains,HashMap<String,Integer> geneNumbers,File out,int draws){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			for(int i=0;i<strains.size();i++){
				Stats s=drawGeneNumbers(strains,geneNumbers,draws,i+1);
				bw.write(s.getAverage()+"\t"+s.getStandardDeviation()+"\n");
				System.out.println("Average of "+(i+1)+" genomes: "+s.getAverage()+" +- "+s.getStandardDeviation());
			}
			bw.close();
		}catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	public static ArrayList<Integer> makeArrayList(int size){
		ArrayList<Integer> temp=new ArrayList<Integer>();
		for(int i=0;i<size;i++){
			temp.add(i);
			
		}
		return temp;
	}
	
	public static Stats drawGroupNumbers(ArrayList<String> strains,HashMap<String,HashMap<Integer,Boolean>> groupNumbers,int draws,int genomes){
		ArrayList<Double> values=new ArrayList<Double>();
		for(int i=0;i<draws;i++){
			HashMap<Integer,Boolean> groupHash=new HashMap<Integer, Boolean>();
			ArrayList<Integer> temp=makeArrayList(strains.size());
			
			for(int j=0;j<genomes;j++){
				int rand=(int)(Math.random()*temp.size());
				groupHash.putAll(groupNumbers.get(strains.get(temp.get(rand))));
				temp.remove(rand);
			}
			values.add(groupHash.size()+0.0);
		}
		return new Stats(values);
	}
	
	public static Stats drawGeneNumbers(ArrayList<String> strains,HashMap<String,Integer> geneNumbers,int draws,int genomes){
		ArrayList<Double> values=new ArrayList<Double>();
		for(int i=0;i<draws;i++){
			double value=0;
			ArrayList<Integer> temp=makeArrayList(strains.size());
			for(int j=0;j<genomes;j++){
				int rand=(int)(Math.random()*temp.size());
				value+=geneNumbers.get(strains.get(temp.get(rand)));
				temp.remove(rand);
			}
			values.add(value);
		}
		return new Stats(values);
	}
	
	public static String getStrainName(String gene){
		StringBuffer strain=new StringBuffer();
		String split2[]=gene.split("-");
		strain.append(split2[0]);
		for(int j=1;j<split2.length-1;j++){
			strain.append("-"+split2[j]);
		}
		return strain.toString();
	}
	
	public static void getGroups(File in,ArrayList<String> strains,HashMap<String,HashMap<Integer,Boolean>> groups){
		try{
			BufferedReader br=new BufferedReader(new FileReader(in));
			String line="";
			while((line=br.readLine())!=null){
				String[] split=line.split("\t");
				int group=Integer.parseInt(split[0]);
				for(int i=1;i<split.length;i++){
					String strain=getStrainName(split[i]);
						if(groups.containsKey(strain)){
							groups.get(strain).put(group,true);
						}else{
							HashMap<Integer,Boolean> temp=new HashMap<Integer, Boolean>();
							temp.put(group,true);
							groups.put(strain,temp);
							strains.add(strain);
							//System.out.println(strain);
						}
				}
			}
			
			br.close();
		}catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	public static void getGeneNumbers(File groups,ArrayList<String> strains,HashMap<String,Integer> numbers){
		try{
			BufferedReader br=new BufferedReader(new FileReader(groups));
			String line="";
			while((line=br.readLine())!=null){
				String[] split=line.split("\t");
				for(int i=1;i<split.length;i++){
					String strain=getStrainName(split[i]);
					if(numbers.containsKey(strain)){
						numbers.put(strain, numbers.get(strain)+1);
					}else{
						numbers.put(strain, 1);
						strains.add(strain);
						//System.out.println(strain);
					}
				}
			}
			
			br.close();
		}catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
	}
}
