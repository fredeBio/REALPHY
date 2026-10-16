package blastTools;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;


public class ReadBlast {
	ArrayList<Integer> startS=new ArrayList<Integer>();
	ArrayList<Integer> endS=new ArrayList<Integer>();
	ArrayList<String> query=new ArrayList<String>();
	ArrayList<String> database=new ArrayList<String>();
	ArrayList<Double> evalue=new ArrayList<Double>();
	ArrayList<Integer> length=new ArrayList<Integer>();
	ArrayList<Integer> startQ=new ArrayList<Integer>();
	ArrayList<Integer> endQ=new ArrayList<Integer>();
	ArrayList<Integer> identities=new ArrayList<Integer>();
	ArrayList<Double> score=new ArrayList<Double>();
	ArrayList<Double> percentIdentity=new ArrayList<Double>();
	
	private void getLines(File f){
		try{	
			BufferedReader br = new BufferedReader(new FileReader(f));
			String line="";
			while((line=br.readLine())!=null){
				String[] split=line.split("\\s+");
				score.add(Double.parseDouble(split[11]));
				startS.add(Integer.parseInt(split[8]));
				endS.add(Integer.parseInt(split[9]));	
				startQ.add(Integer.parseInt(split[6]));
				endQ.add(Integer.parseInt(split[7]));
				query.add(split[0]);
				evalue.add(Double.parseDouble(split[10]));
				database.add(split[1]);
				length.add(Integer.parseInt(split[3]));
				Double d=(Double.parseDouble(split[2])/100)*length.get(length.size()-1);
				identities.add(d.intValue());
				percentIdentity.add(Double.parseDouble(split[2]));
			}
			br.close();
		}catch(IOException e){
			System.err.println(e.toString());
		}
	}
	
	public ReadBlast(File f){
		getLines(f);
	}
	public ArrayList<String> getQuery(){
		return query;
	}
	
	public ArrayList<Double> getEvalue(){
		return evalue;
	}
	
	public ArrayList<String> getDatabase(){
		return database;
	}
	public ArrayList<Integer> getIdentities(){
		return identities;
	}
	public ArrayList<Integer> getLength(){
		return length;
	}
	public ArrayList<Integer> getStartQuery(){
		return startQ;
	}
	
	public ArrayList<Integer> getEndQuery(){
		return endQ;
	}
	public ArrayList<Integer> getStartDB(){
		return startS;
	}
	public ArrayList<Double> getScore(){
		return score;
	}
	public ArrayList<Double> getPercentIdentity(){
		return percentIdentity;
	}
	public ArrayList<Integer> getEndDB(){
		return endS;
	}
	public String get(int i){
		return query.get(i)+"\t"+database.get(i)+"\t"+identities.get(i)+"\t"+startQ.get(i)+"\t"+endQ.get(i)+"\t"+startS.get(i)+"\t"+endS.get(i)+"\t"+evalue.get(i)+"\t"+score.get(i);
	}
}
