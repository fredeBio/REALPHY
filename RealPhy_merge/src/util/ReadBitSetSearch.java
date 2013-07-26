package util;

import java.util.ArrayList;


public class ReadBitSetSearch implements ReadSearchOutput{
	ArrayList<Integer> end=new ArrayList<Integer>();
	ArrayList<Integer> start=new ArrayList<Integer>();
	ArrayList<String> id=new ArrayList<String>();
	
	public ReadBitSetSearch(ArrayList<Integer> pos,int size){
		generateIDs(pos);
		generateEnds(pos,size);
		generateStarts(pos,size);
	}
	private  void generateIDs(ArrayList<Integer> pos) {
		System.err.println("ALL positions have the same query ID!");
		for(int i=0;i<pos.size();i++){
			id.add("0");
		}
		
	}
	private  void generateStarts(ArrayList<Integer> pos,int size) {
		
		for(int i=0;i<pos.size();i++){
			if(pos.get(i)>0)start.add(pos.get(i));
			else start.add(Math.abs(pos.get(i))+size);
		}
		
	}
	private  void generateEnds(ArrayList<Integer> pos,int size) {
		
		for(int i=0;i<pos.size();i++){
			if(pos.get(i)>0)end.add(pos.get(i)+size);
			else end.add(Math.abs(pos.get(i)));
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
		return id;
	}

	@Override
	public ArrayList<Integer> getStart() {
		// TODO Auto-generated method stub
		return start;
	}

}
