package util;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map.Entry;

public class HashMapManip extends HashMap<Integer, Integer>{
	static final long serialVersionUID=1;
	public  void add(HashMap<Integer,Integer> h){
		
		Iterator<Entry<Integer,Integer>> it=h.entrySet().iterator();
		while(it.hasNext()){
			Entry<Integer,Integer> e=it.next();
			if(containsKey(e.getKey())){
				put(e.getKey(),get(e.getKey())+e.getValue());
			}else{
				put(e.getKey(),e.getValue());
			}
		}
		
		
	}
}
