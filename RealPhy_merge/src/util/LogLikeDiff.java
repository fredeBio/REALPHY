package util;

public class LogLikeDiff {
	double FullCorr;
	double FullInc;
	double RPCorr;
	double RPInc;
	double Total;
	public static final int FC=0;
	public static final int FI=1;
	public static final int RC=2;
	public static final int RI=3;
	public static final int TO=4;
	public LogLikeDiff(double total,double fullCorr,double fullInc,double rpCorr,double rpInc){
		FullCorr=fullCorr;
		FullInc=fullInc;
		RPCorr=rpCorr;
		RPInc=rpInc;
		Total=total;
	}
	public double get(int variable){
		if(variable==FC){
			return FullCorr;
		}
		if(variable==FI){
			return FullInc;
		}
		if(variable==RC){
			return RPCorr;
		}
		if(variable==RI){
			return RPInc;
		}
		if(variable==TO){
			return Total;
		}
		return Double.NaN;
	}
}
