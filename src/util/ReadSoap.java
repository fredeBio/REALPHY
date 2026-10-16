package util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;


public class ReadSoap implements ReadSearchOutput{
	ArrayList<Integer> start;
	ArrayList<Integer> end;
	ArrayList<String> query;
	ArrayList<String> database;
	ArrayList<Boolean> forward;
	public ReadSoap(File f){
		start=new ArrayList<Integer>();
		end=new ArrayList<Integer>();
		query=new ArrayList<String>();
		database=new ArrayList<String>();
		forward=new ArrayList<Boolean>();
		readLines(f);
	}
	private void readLines(File f){
		try{
			BufferedReader br=new BufferedReader(new FileReader(f));
			String line="";
			while((line=br.readLine())!=null){
				String splits[]=line.split("\\s+");
				start.add(Integer.parseInt(splits[8]));
				end.add(Integer.parseInt(splits[8])+Integer.parseInt(splits[5])-1);
				query.add(splits[0]);
				database.add(splits[7]);
				forward.add(splits[6].equals("+"));
			}
		}catch(IOException e){
			System.err.println(e.toString());
		}
	}
	
	public ArrayList<String> getQuery(){
		return query;
	}
	
	
	public ArrayList<String> getDatabase(){
		return database;
	}
	
	public ArrayList<Integer> getStart(){
		return start;
	}
	
	public ArrayList<Integer> getEnd(){
		return end;
	}
	public ArrayList<Boolean> getForward(){
		return forward;
	}
}
