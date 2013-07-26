package pairwiseAlignment;

public class Alignment {
	String alA;
	String alB;
	public Alignment(String a,String b){
		alA=a;
		alB=b;
	}
	public String toString(){
		return alA+"\n"+alB+"\n";
	}
}
