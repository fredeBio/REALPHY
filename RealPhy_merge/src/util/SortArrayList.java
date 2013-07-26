package util;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.TreeMap;
import java.util.Map.Entry;

public  class SortArrayList {

	
	public static <T> ArrayList<T> sort(ArrayList<T> list){
		TreeMap<T,Integer> tm=new TreeMap<T, Integer>();
		ArrayList<T> sorted=new ArrayList<T>();
		for(int i=0;i<list.size();i++){
			if(!tm.containsKey(list.get(i))){
				tm.put(list.get(i),1);
			}else{
				tm.put(list.get(i), tm.get(list.get(i))+1);
			}
		}
		Iterator<Entry<T,Integer>> it=tm.entrySet().iterator();
		while(it.hasNext()){
			Entry<T,Integer> e=it.next();
			int rep=e.getValue();
			T key=e.getKey();
			for(int i=0;i<rep;i++){
				sorted.add(key);
			}
		}
		return sorted;
	}
	
	
	
}
