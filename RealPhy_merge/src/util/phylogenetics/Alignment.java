package util.phylogenetics;

import java.util.*;

import util.*;

public class Alignment {
	int length;
	int numSeqs;
	ArrayList<StringBuffer> seqs=new ArrayList<StringBuffer>();
	ArrayList<String> idents=new ArrayList<String>();
	public ArrayList<String> getIdents(){
		return idents;
	}
	public Alignment(ArrayList<Fasta> fas)throws RuntimeException{
		numSeqs=fas.size();
		idents=new ArrayList<String>();
		seqs=new ArrayList<StringBuffer>();
		if(numSeqs>1)checkLengthAndBases(fas);
		else throw new RuntimeException("The alignment needs to contain at least two sequences.\n");
	}
	
	public Alignment(Alignment al){
		numSeqs=al.numSeqs;
		length=al.length;
		for(int i=0;i<al.numSeqs;i++){
			seqs.add(new StringBuffer(al.seqs.get(i)));
		}
		for(int i=0;i<al.numSeqs;i++){
			idents.add(new String(al.idents.get(i)));
		}
	}
	public Alignment(Alignment al,int start ,int end){
		numSeqs=al.numSeqs;
		length=al.length;
		for(int i=0;i<al.numSeqs;i++){
			seqs.add(new StringBuffer(al.seqs.get(i).substring(start, end)));
		}
		for(int i=0;i<al.numSeqs;i++){
			idents.add(new String(al.idents.get(i)));
		}
	}

	public Alignment(){
		numSeqs=0;
		length=0;
		seqs=new ArrayList<StringBuffer>();
		idents=new ArrayList<String>();
	}
	
	public ArrayList<Fasta> toFasta(){
		ArrayList<Fasta> fas=new ArrayList<Fasta>();
		for(int i=0;i<numSeqs;i++){
			fas.add(new Fasta(idents.get(i),seqs.get(i).toString()));
		}
		return fas;
	}
	
	public void addColumn(String col) {
		if(numSeqs==0){
			numSeqs=col.length();
			length=0;
			for(int i=0;i<numSeqs;i++){
				idents.add(i+"");
				seqs.add(new StringBuffer());
			}
		}else{
			if(col.length()!=numSeqs){
				throw new RuntimeException("Column length ("+col.length()+") has to be identical to the existing sequence number ("+numSeqs+") in order to be added.");
			}
		}
		
		col=col.toUpperCase();
		for(int i=0;i<numSeqs;i++){
			seqs.get(i).append(checkBase(col.charAt(i),length,i));
		}
		length++;
	}
	
	public void changeIdents(ArrayList<String> Idents){
		idents=Idents;
	}
	
	public char getBase(int seq,int pos){
		return seqs.get(seq).charAt(pos);
	}
	public String getSequence(int num){
		return seqs.get(num).toString();
	}
	public String getIdent(int num){
		return idents.get(num);
	}
	public int getLength(){
		return length;
	}
	public int getNumSeqs(){
		return numSeqs;
	}
	
	public String getColumn(int pos){
		StringBuffer sb=new StringBuffer();
		for(int i=0;i<numSeqs;i++){
			sb.append(seqs.get(i).charAt(pos));
		}
		return sb.toString();
	}
	
	public Alignment getRegion(int start,int end){
		Alignment al=new Alignment(this,start,end);

		return al;
	}
	
	private void checkLengthAndBases(ArrayList<Fasta> fas)throws RuntimeException{
		length=fas.get(0).getSequence().length();
		for(int i=0;i<numSeqs;i++){
			if(fas.get(i).getSequence().length()!=length){
				throw new RuntimeException("All sequences have to be of the same length.\n"+fas.get(0).getIdent()+" length: "+length+"\n"+fas.get(i).getIdent()+" length: "+fas.get(i).getSequence().length()+"\n");
			}
			StringBuffer seq=new StringBuffer(fas.get(i).getSequence().toUpperCase());
			idents.add(fas.get(i).getIdent());
			for(int j=0;j<length;j++){
				char c=seq.charAt(j);
				seq.setCharAt(j, checkBase(c,j,i));
			}
			seqs.add(seq);

		}
	}
	private char checkBase(char c,int pos,int seqNumber){
		if(c!='A'&&c!='T'&&c!='C'&&c!='G'&&c!='-'){
			System.err.println("Changed character "+c+" in sequence "+idents.get(seqNumber)+" at position "+pos+" to "+"-.");

			return '-';
		}else{
			return c;
		}
	}
}
