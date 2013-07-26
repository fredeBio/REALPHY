package util;

import java.util.Map.Entry;

public class EntryExpanded<K,V> implements Entry<K, V> {

	private K Key;
	private V Value;
	@Override
	public K getKey() {
		// TODO Auto-generated method stub
		return Key;
	}
	@Override
	public V getValue() {
		// TODO Auto-generated method stub
		return Value;
	}
	@Override
	public V setValue(V value) {
		// TODO Auto-generated method stub
		Value=value;
		return Value;
	}
	public EntryExpanded(K key,V value){
		Key=key;
		Value=value;
	}
	
}
