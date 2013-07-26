package util;

import java.util.ArrayList;

import statistics.Stats;

public class BinLength {
	ArrayList<ArrayList<Double>> bins=new ArrayList<ArrayList<Double>>();
	boolean binned=false;
	boolean binStats=false;
	ArrayList<Double> binContentHistX=new ArrayList<Double>();
	ArrayList<Double> binContentHistY=new ArrayList<Double>();
	Stats binContentHistStats;
	double binSize;
	int NumBins;
	ArrayList<Double> list;
	
	public BinLength(ArrayList<Double> List){
		list=List;
	}
	private void initBin(ArrayList<ArrayList<Double>> bins){
		for(int i=0;i<NumBins;i++){
			bins.add(new ArrayList<Double>());
		}
		
	}
	private void calculateBinsStats(){
		binStats=true;
		int size=bins.size();
		for(int i=0;i<size;i++){
			int numElements=bins.get(i).size();
			binContentHistY.add(numElements+0.0);
			binContentHistStats=new Stats(binContentHistY);
		}
	}
	
	private void makeXValues(double min,double binSize){
		for(int i=0;i<NumBins;i++){
			binContentHistX.add(min+i*binSize);
		}
	}
	private void binData(){
			binned=true;
			double min=Stats.getMin(list);
			double max=Stats.getMax(list);
			
			binSize=(max-min)/NumBins;
			initBin(bins);
			makeXValues(min,binSize);
			for(int i=0;i<list.size();i++){
				double entry=list.get(i);
				int bin=(int)((entry-min)/binSize);
				if(bin==NumBins){
					bin=bin-1;
				}
				bins.get(bin).add(list.get(i));
			}
			calculateBinsStats();
	}
	public Stats getBinContentStatsLength(int numBins){
		if(binStats&&numBins==NumBins){
			return binContentHistStats;
		}else{
			NumBins=numBins;
			binData();
			return binContentHistStats;
		}
		
	}
	public ArrayList<Double> getBinsSizesPositionX(int numBins){
		if(!binned)runBinning(numBins);
		return binContentHistX;
	}
	
	public ArrayList<Double> getBinSizes(int numBins){
		if(!binned)runBinning(numBins);
		return binContentHistY;
	}
	private void runBinning(int numBins){
		if(!binned||(numBins!=NumBins)){
			NumBins=numBins;
			binData();
		}
	}
}
