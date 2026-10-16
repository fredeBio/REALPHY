package util;

import java.io.Serializable;

public  class Interval<T> implements Comparable<Interval<T>>, Serializable{
	Integer start;
	Integer end;
	public Interval(int s,int e){
		start=s;
		end=e;
	}
	
	public int getStart(){
		return start;
	}
	
	public int getEnd(){
		return end;
	}
	
	public int compareTo(Interval<T> o1) {
		// TODO Auto-generated method stub
		if(start.compareTo(o1.start)!=0){
			return start.compareTo(o1.start);
		}else 
			return end.compareTo(o1.end);

	}
	public boolean overlapsWith(Interval<T> other){
		if(start<=other.end && start>=other.start){
			return true;
		}else if(end<=other.end && end>=other.start){
			return true;
		}else if(end>=other.end && start<=other.end){
			return true;
		}else if(end>=other.start && start<=other.start){
			return true;
		}
		return false;
	}
	public  void append(T interval){
		
	}
}
