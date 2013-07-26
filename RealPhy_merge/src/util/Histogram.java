package util;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.TreeMap;
import java.util.Map.Entry;
public class Histogram <T>{
	TreeMap<T,Double> histogram=new TreeMap<T, Double>();
	ArrayList<TreeMap<T,Double>> elementHash=new ArrayList<TreeMap<T,Double>>();
	TreeMap<T,Double> avg=new TreeMap<T, Double>();
	TreeMap<T,Double> stdv=new TreeMap<T, Double>();

	int numHist=0;
	public Histogram(T[] elements){
		makeHistogram(elements);
		numHist=1;
	}

	public Histogram(TreeMap<T,Double> h){
		histogram=h;
	}
	
	//divides all the elements in the histogram by the number of histograms added to achieve an average value
	public void setAverage(){
		Iterator<Entry<T,Double>> it=histogram.entrySet().iterator();
		while(it.hasNext()){
			Entry<T,Double> e=it.next();
			avg.put(e.getKey(),e.getValue()/numHist);
			
		}
	}
	
	public void setStandardDeviation(){
		HashMap<T,Double> help=new HashMap<T, Double>();
		for(int i=0;i<elementHash.size();i++){
			Iterator<Entry<T,Double>> it=elementHash.get(i).entrySet().iterator();
			while(it.hasNext()){
				Entry<T,Double> e=it.next();
				if(!help.containsKey(e.getKey())){
					help.put(e.getKey(), Math.pow(e.getValue()-avg.get(e.getKey()),2)/elementHash.size());
				}else{
					help.put(e.getKey(), help.get(e.getKey())+Math.pow(e.getValue()-avg.get(e.getKey()),2)/elementHash.size());
				}
			}
		}
		Iterator<Entry<T,Double>> it=help.entrySet().iterator();
		while(it.hasNext()){
			Entry<T,Double> e=it.next();
			stdv.put(e.getKey(),Math.sqrt(e.getValue()));
			
		}
	}
	
	public void add(T[] elements){

		elementHash.add(makeHistogram(elements));
		numHist++;
	}
	public void max(T[] elements){
		HashMap<T,Double> H2=new HashMap<T, Double>();
		for(int i=0;i<elements.length;i++){
			if(H2.containsKey(elements[i])){
				H2.put(elements[i],H2.get(elements[i])+1);
			}else{
				H2.put(elements[i], 1.0);
			}
		}
		Iterator<Entry<T,Double>> it=H2.entrySet().iterator();
		while(it.hasNext()){
			Entry<T,Double> e=it.next();
			if(histogram.containsKey(e.getKey())){
				if(histogram.get(e.getKey())<e.getValue()){
					histogram.put(e.getKey(), e.getValue());
				}
			}else{
				histogram.put(e.getKey(),e.getValue());
			}
		}
	}
	private TreeMap<T,Double> makeHistogram(T[] elements){
		TreeMap<T,Double> newHash=new TreeMap<T, Double>();
		for(int i=0;i<elements.length;i++){
			if(histogram.containsKey(elements[i])){
				histogram.put(elements[i],histogram.get(elements[i])+1);
			}else{
				histogram.put(elements[i], 1.0);
			}
			if(newHash.containsKey(elements[i])){
				newHash.put(elements[i],newHash.get(elements[i])+1);
			}else{
				newHash.put(elements[i], 1.0);
			}
		}
		return newHash;
	}
	
	public TreeMap<T,Double> getHistogram(){
		return histogram;
	}
	
	public Double[] getValues(){
		return histogram.values().toArray(new Double[0]);
	}
	
	public TreeMap<Double,ArrayList<T>> sortByValues(){
		TreeMap<Double, ArrayList<T>> sorted=new TreeMap<Double, ArrayList<T>>();
		Iterator<Entry<T,Double>> it=histogram.entrySet().iterator();
		while(it.hasNext()){
			Entry<T,Double> e=it.next();
			if(sorted.containsKey(e.getValue())){
				sorted.get(e.getValue()).add(e.getKey());
			}else{
				ArrayList<T> a=new ArrayList<T>();
				a.add(e.getKey());
				sorted.put(e.getValue(),a);
			}
		}
		return sorted;
		
	}
	
	public static TreeMap<Integer,Double> calculateProportionIntegral(TreeMap<Integer,Double> histogram){
		Iterator<Entry<Integer,Double>> it=histogram.entrySet().iterator();
		TreeMap<Integer,Double> tm=new TreeMap<Integer, Double>();
		double sum=0;
		double prop=0;
		while(it.hasNext()){
			Entry<Integer,Double> e=it.next();
			sum+=e.getKey()*e.getValue();
		}
		it=histogram.entrySet().iterator();
		while(it.hasNext()){
			Entry<Integer,Double> e=it.next();
			prop+=e.getKey()*e.getValue();
			tm.put(e.getKey(), prop/sum);
		}
		return tm;
	}
	
	public void set(TreeMap<T,Double> h){
		histogram=h;
	}
	
	public void writeSortedByValues(File out){
		TreeMap<Double, ArrayList<T>> sorted=sortByValues();
		Iterator<Entry<Double,ArrayList<T>>> it=sorted.descendingMap().entrySet().iterator();
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			while(it.hasNext()){
				Entry<Double,ArrayList<T>> e=it.next();
				for(int j=0;j<e.getValue().size();j++){
					bw.write(e.getValue().get(j)+"\t"+e.getKey()+"\n");
				}

			}
			bw.close();
		}catch(IOException e){
			System.err.println(e.toString());
		}
		
		
		
	}
	
	public TreeMap<T,Double> getStandardDeviation(){
		return stdv;
	}
	
	public TreeMap<T,Double> getAverage(){
		return avg;
	}
	
	public String write(String what){
		StringBuffer out=new StringBuffer();
		Iterator<Entry<T,Double>> it=what.equalsIgnoreCase("std")?histogram.entrySet().iterator():what.equalsIgnoreCase("stdv")?stdv.entrySet().iterator():what.equalsIgnoreCase("avg")?avg.entrySet().iterator():null;
		int old=0;
		while(it.hasNext()){
			Entry<T,Double> e=it.next();
			if(e.getKey() instanceof Integer && (Integer)e.getKey()<500){
				for(int i=1;i<(Integer)e.getKey()-old;i++)out.append((old+i)+"\t"+0+"\n");
				old=(Integer)e.getKey();
			}
			out.append(e.getKey()+"\t"+e.getValue()+"\n");
		}
		return out.toString();
	}
	
	public void write(File out,String what){
		
		try{
			BufferedWriter bw = new BufferedWriter(new FileWriter(out));
			bw.write(write(what));
			bw.close();
		}catch(IOException e){
			System.err.println(e.toString());
		}
	
	}
	public static String write(HashMap<Integer,Integer> hist,int accuracy){
		StringBuffer out=new StringBuffer();
		int i=0;
		while(hist.size()>0 && i<5000000){
			if(hist.containsKey(i)){
				out.append(i+"\t"+hist.get(i)+"\n");
				hist.remove(i);
			}
			i+=accuracy;
			if(i % 10000000 == 0){
				System.out.println(i);
			}
		}
		return out.toString();
	}
	//what is either stdv for standard deviation, std for accumulated data values (either through max or add), or avg for average over all datavalues
	public static String write(HashMap<String,Integer> hist){
		StringBuffer out=new StringBuffer();
		Iterator<Entry<String,Integer>> it=hist.entrySet().iterator();
		while(it.hasNext()){
			Entry<String,Integer> e=it.next();
			out.append(e.getKey()+"\t"+e.getValue()+"\n");
		}
			
		return out.toString();
	}
	public static void write(HashMap<String,Integer> hist,File f){
		try{
			BufferedWriter bw = new BufferedWriter(new FileWriter(f));
			bw.write(write(hist));
			bw.close();
		}catch(IOException e){
			System.err.println(e.toString());
		}
	}
	
	public static void write(HashMap<Integer,Integer> hist,int accuracy,File f){
		try{
			BufferedWriter bw = new BufferedWriter(new FileWriter(f));
			bw.write(write(hist,accuracy));
			bw.close();
		}catch(IOException e){
			System.err.println(e.toString());
		}
	}
	

	
}
