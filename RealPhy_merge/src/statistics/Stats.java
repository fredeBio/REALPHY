package statistics;

import java.util.*;

public class Stats {
	
	public class CoupledEntry implements Comparable<CoupledEntry>{
		double x;
		double y;
		public int compareTo(CoupledEntry e){
			return x<e.x?-1:x==e.x?0:1;
		}
		public CoupledEntry(double X,double Y){
			x=X;
			y=Y;
		}
	}
	
	private ArrayList<Double> listX=new ArrayList<Double>();
	private ArrayList<Double> listY=new ArrayList<Double>();
	ArrayList<ArrayList<Double>> binsX=new ArrayList<ArrayList<Double>>();
	ArrayList<ArrayList<Double>> binsY=new ArrayList<ArrayList<Double>>();

	ArrayList<Double> binsYMean=new ArrayList<Double>();
	ArrayList<Double> binsYStdv=new ArrayList<Double>();
	ArrayList<Double> binsYStdErr=new ArrayList<Double>();
	ArrayList<Double> binsXMean=new ArrayList<Double>();
	private double stdDev;
	private double avg;
	private boolean meanCalc;
	private boolean stdvCalc;
	boolean twoLists=false;
	boolean binned=false;

	double binSize;
	int NumBins;

	
	public Stats(ArrayList<Double> List){
		listX=List;
		meanCalc=false;
		stdvCalc=false;

		
	}
	public ArrayList<Double> getList(){
		return listX;
	}
	
	public static ArrayList<Double> toArrayList(double[] list){
		ArrayList<Double> newList=new ArrayList<Double>();
		for(int i=0;i<list.length;i++){
			newList.add(list[i]);
		}
		return newList;
	}
	
	public Stats(double[] List){
		listX=toArrayList(List);
		meanCalc=false;
		stdvCalc=false;

		
	}
	
	public Stats(ArrayList<Double> x,ArrayList<Double> y){
		listX=x;
		listY=y;
		meanCalc=false;
		stdvCalc=false;
		twoLists=true;
		
	}
    public static double getMin(double[] all){
    	double min=all[0];
    	for(int i=0;i<all.length;i++){
    		if(min>all[i]){
    			min=all[i];
    		}
    	}
    	return min;
    }
    public static double getMin(ArrayList<Double> all){
    	double min=all.get(0);
    	for(int i=0;i<all.size();i++){
    		if(min>all.get(i)){
    			min=all.get(i);
    		}
    	}
    	return min;
    }
    public static double getMax(double[] all){
    	double max=all[0];
    	for(int i=0;i<all.length;i++){
    		if(max<all[i]){
    			max=all[i];
    		}
    	}
    	return max;
    }
    public static double getMax(ArrayList<Double> all){
    	double max=all.get(0);
    	for(int i=0;i<all.size();i++){
    		if(max<all.get(i)){
    			max=all.get(i);
    		}
    	}
    	return max;
    }
	
    private ArrayList<CoupledEntry> makeList(){
    	ArrayList<CoupledEntry> ce=new ArrayList<Stats.CoupledEntry>();
    	for(int i=0;i<listX.size();i++){
    		ce.add(new CoupledEntry(listX.get(i),listY.get(i)));
    	}
    	return ce;
    }
    

    
    private void sortLists(){
    	ArrayList<CoupledEntry> al=makeList();
    	Collections.sort(al);
    	listX=new ArrayList<Double>();
    	listY=new ArrayList<Double>();

    	for(int i=0;i<al.size();i++){
    		listX.add(al.get(i).x);
    		listY.add(al.get(i).y);
    	}
    }
    

    

	
	private void binData(){
		if(twoLists){
			binned=true;
//			double min=getMin(listX);
//			double max=getMax(listX);
			binSize=listX.size()/NumBins;
			int j=-1;
			sortLists();
			
			for(int i=0;i<listX.size();i++){
				if(i>=(j+1)*binSize){
					binsX.add(new ArrayList<Double>());
					binsY.add(new ArrayList<Double>());
					j++;
				}
				binsX.get(j).add(listX.get(i));
				binsY.get(j).add(listY.get(i));
			}
			calculateBinsStats();
		}else{
			System.err.println("For a single list binning is not implemented yet.");
		}
	}
	
	
	private void calculateBinsStats(){
		int size=binsY.size();
		for(int i=0;i<size;i++){
			Stats s=new Stats(binsY.get(i));
			binsYMean.add(s.getAverage());
			binsYStdv.add(s.getStandardDeviation());
			binsYStdErr.add(s.getStandardError());
			s=new Stats(binsX.get(i));
			binsXMean.add(s.getAverage());
		}
	}
	

	

	

	

	
	public ArrayList<Double> getBinsYStdDev(int numBins){
		runBinning(numBins);
		return binsYStdv;
	}
	
	public ArrayList<Double> getBinsYStdErr(int numBins){
		runBinning(numBins);
		return binsYStdErr;
	}
	
	public ArrayList<Double> getBinsYMean(int numBins){
		runBinning(numBins);
		return binsYMean;
	}
	
	public ArrayList<Double> getBinsXMean(int numBins){
		runBinning(numBins);
		return binsXMean;
	}
	
	public ArrayList<ArrayList<Double>> getBinsX(int numBins){
		runBinning(numBins);
		return binsX;
	}

	public ArrayList<ArrayList<Double>> getBinsY(int numBins){
		runBinning(numBins);

		return binsY;
	}
	
	public ArrayList<ArrayList<Double>> getBinsX(){
		
		return binsX;
	}

	public ArrayList<ArrayList<Double>> getBinsY(){
		

		return binsY;
	}
	
	private void runBinning(int numBins){
		if(!binned&&twoLists){
			NumBins=numBins;
			binData();
		}else if(binned&&twoLists&&numBins!=NumBins){
			NumBins=numBins;
			binData();
		}
	}
	


	
	public double getBinSize(){
		if(binned){
			return binSize;
		}else{
			System.err.println("Data has not been binned yet. Please run binData(double numBins) first.");
			return -1;
		}
	}
	public ArrayList<Integer> getIndexListBelowDoubleStdv(){
		ArrayList<Integer> indeces=new ArrayList<Integer>();
		if(!meanCalc){
			calculateMean();
			meanCalc=true;
		}
		if(!stdvCalc){
			calculateStdDev();
			stdvCalc=true;
		}
		for(int i=0;i<listX.size();i++){
			if(listX.get(i)<avg-stdDev*2){
				indeces.add(i);
			}
		}
		return indeces;
	}
	public ArrayList<Integer> getIndexListAboveDoubleStdv(){
		ArrayList<Integer> indeces=new ArrayList<Integer>();
		if(!meanCalc){
			calculateMean();
			meanCalc=true;
		}
		if(!stdvCalc){
			calculateStdDev();
			stdvCalc=true;
		}
		for(int i=0;i<listX.size();i++){
			if(listX.get(i)>stdDev*2+avg){
				indeces.add(i);
			}
		}
		return indeces;
	}
	
	public double unexplainedVar(){
		double avg=getAverage();
		double sum=0;
		for(int i=0;i<listX.size();i++){
			sum+=Math.pow(avg-listX.get(i),2);
		}
		return sum/(listX.size()-2);
	}
	
	private void calculateMean(){
		double sum=0;
		for(int i=0;i<listX.size();i++){
			sum+=listX.get(i);
		}
		avg=sum/listX.size();
	}
	
	private void calculateStdDev(){
		double help=0;
		for(int i=0;i<listX.size();i++){
			help+=Math.pow((listX.get(i)-avg),2);
		}
		stdDev=Math.sqrt(help/listX.size());
	}
	public double getStandardError(){
		if(!meanCalc){
			calculateMean();
			calculateStdDev();
			meanCalc=true;
			stdvCalc=true;
		}else if(!stdvCalc){
			calculateStdDev();
			stdvCalc=true;
		}
		return stdDev/Math.sqrt(listX.size());
	}
	public double getStandardDeviation(){
		if(!meanCalc){
			calculateMean();
			meanCalc=true;
		}
		if(!stdvCalc){
			calculateStdDev();
			stdvCalc=true;
		}
		return stdDev;
	}
	
	public double getAverage(){
		if(!meanCalc){
			meanCalc=true;
			calculateMean();	
		}
		return avg;
	}
}
