package util;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.TreeMap;
public class SortPositions  implements Comparator<Integer>{
	ArrayList<Integer> sorted;
	public SortPositions(ArrayList<Integer> list){
		sorted=new ArrayList<Integer>();
		sortInt(list);
	}
	
	public ArrayList<Integer> getList(){
		return sorted;
	}
	
	private void sortInt(ArrayList<Integer> list){
		TreeMap<Integer,Boolean> tm=new TreeMap<Integer, Boolean>(this);
		for(int i=0;i<list.size();i++){
			tm.put(list.get(i),true);
		}
		Iterator<Integer> it=tm.keySet().iterator();
		while(it.hasNext()){
			sorted.add(it.next());
		}
	}

	public static ArrayList<Integer> sort(ArrayList<Integer> list){
		TreeMap<Integer,Boolean> tm=new TreeMap<Integer, Boolean>();
		ArrayList<Integer> sorted=new ArrayList<Integer>();
		for(int i=0;i<list.size();i++){
			tm.put(list.get(i),true);
		}
		Iterator<Integer> it=tm.keySet().iterator();
		while(it.hasNext()){
			sorted.add(it.next());
		}
		return sorted;
	}
	
	@Override
	public int compare(Integer o1, Integer o2) {
		if(Math.abs(o1)>Math.abs(o2)){
			return 1;
		}else if(Math.abs(o2)>Math.abs(o1)){
			return -1;
		}else
			return 0;
	}
	
}
