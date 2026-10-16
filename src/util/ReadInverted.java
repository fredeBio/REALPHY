package util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class ReadInverted implements ReadSearchOutput{
	ArrayList<Integer> end=new ArrayList<Integer>();
	ArrayList<Integer> start=new ArrayList<Integer>();
	ArrayList<String> query=new ArrayList<String>();
	public ReadInverted(File f){
		read(f);
	}
	
	private void read(File f){
		try{
			BufferedReader br=new BufferedReader(new FileReader(f));
			String line="";
			int i=0;
			while((line=br.readLine())!=null){
				i++;
				if(!line.matches("^\\d.+"))continue;
				String[] split=line.split("\\s+");
				start.add(Integer.parseInt(split[0]));
				end.add(Integer.parseInt(split[1]));
				query.add("IV"+i);
			}
			
			
		}catch(IOException e){
			System.err.println(e.toString());
		}
		
		
	}


	@Override
	public ArrayList<Integer> getEnd() {
		// TODO Auto-generated method stub
		return end;
	}

	@Override
	public ArrayList<String> getQuery() {
		// TODO Auto-generated method stub
		return query;
	}

	@Override
	public ArrayList<Integer> getStart() {
		// TODO Auto-generated method stub
		return start;
	}	
}