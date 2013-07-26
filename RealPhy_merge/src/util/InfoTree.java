package util;

import java.util.ArrayList;
import java.util.HashMap;


public class InfoTree extends IntervalTree<Info> {
    
	 public void parseTree(IntervalNode<Info> t,ArrayList<Info> intervals,HashMap<String,Boolean> feature){
	    	if(t==nullNode)return; 
	    	parseTree(t.left,intervals,feature);
	       	if(feature.containsKey(t.element.getFeature())){
	    		intervals.add(new Info(t.element));
	    	}
	    	parseTree(t.right,intervals,feature);
	    }
	    public void parseTree(IntervalNode<Info> t,ArrayList<Info> intervals,HashMap<String,Boolean> feature,int minLength){
	    	if(t==nullNode)return; 
	    	parseTree(t.left,intervals,feature,minLength);
	       	if(feature.containsKey(t.element.getFeature()) && (t.element.end-t.element.start)>=minLength){
	    		intervals.add(new Info(t.element));
	    	}
	    	parseTree(t.right,intervals,feature,minLength);
	    }
	
	public void parseTree(IntervalNode<Info> t,ArrayList<Info> intervals,String feature){
    	if(t==nullNode)return; 
    	parseTree(t.left,intervals,feature);
       	if(t.element.getFeature().equalsIgnoreCase(feature)){
    		intervals.add(new Info(t.element));
    	}
    	parseTree(t.right,intervals,feature);
    }
    public void parseTree(IntervalNode<Info> t,ArrayList<Info> intervals,String feature,int minLength){
    	if(t==nullNode)return; 
    	parseTree(t.left,intervals,feature,minLength);
       	if(t.element.getFeature().equalsIgnoreCase(feature) && (t.element.end-t.element.start)>=minLength){
    		intervals.add(new Info(t.element));
    	}
    	parseTree(t.right,intervals,feature,minLength);
    }
    public void parseTree(IntervalNode<Info> t,ArrayList<Info> intervals){
    	if(t==nullNode)return; 
    	parseTree(t.left,intervals);
    	intervals.add(new Info(t.element));    
    	parseTree(t.right,intervals);
    }
    
    public boolean checkOverlap(Info inf){
    	ArrayList<Info> al=new ArrayList<Info>();
    	this.search(inf, al);
    	if(al.size()>0){
    		return true;
    	}
    	return false;
    }
    
    public ArrayList<Info> getOverlap(Info inf){
    	ArrayList<Info> al=new ArrayList<Info>();
    	this.search(inf, al);
    	
    	return al;
    }
    
    public int getOverlapLength(Info inf){
    	ArrayList<Info> al=new ArrayList<Info>();
    	this.search(inf, al);
    	int[] array=new int[inf.size()];
    	int infStart=inf.getStart();
    	int infEnd=inf.getEnd();
    	for(int i=0;i<inf.size();i++){
    		array[i]=0;
    	}
    	for(int i=0;i<al.size();i++){
    		Info temp=al.get(i);
    		for(int j=temp.getStart();j<temp.getEnd();j++){
    			if(j>=infStart&&j<=infEnd){
    				array[j-infStart]=1;
    			}
    		}
    	}
    	int sum=0;
    	for(int i=0;i<array.length;i++){
    		if(array[i]==1){
    			sum++;
    		}
    	}
    	
    	return sum;
    }
    
    public ArrayList<Info> parseTree(){
    	ArrayList<Info> intervals=new ArrayList<Info>();
    	parseTree(root,intervals);
    	return intervals;
    }
    public ArrayList<Info> parseTree(String feature){
    	ArrayList<Info> intervals=new ArrayList<Info>();
    	parseTree(root,intervals,feature);
    	return intervals;
    }
    public ArrayList<Info> parseTree(HashMap<String,Boolean> feature){
    	ArrayList<Info> intervals=new ArrayList<Info>();
    	parseTree(root,intervals,feature);
    	return intervals;
    }
    public ArrayList<Info> parseTree(String feature,int minLength){
    	ArrayList<Info> intervals=new ArrayList<Info>();
    	parseTree(root,intervals,feature,minLength);
    	return intervals;
    }
    public ArrayList<Info> parseTree(HashMap<String,Boolean> feature,int minLength){
    	ArrayList<Info> intervals=new ArrayList<Info>();
    	parseTree(root,intervals,feature,minLength);
    	return intervals;
    }
	public ArrayList<Info> delete(ArrayList<Info> genesFlanked) {
		for(int i=0;i<genesFlanked.size();i++){
			remove(genesFlanked.get(i));
		}
		
		return parseTree();
	}
	public InfoTree(){
		super();
	}
	public InfoTree(ArrayList<Info> infos){
		for(int i=0;i<infos.size();i++){
			insert(infos.get(i));
		}
	}
	
	
}
