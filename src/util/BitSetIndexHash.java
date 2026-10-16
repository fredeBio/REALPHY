package util;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
//orientation = true -> reverse positions are negative numbers
public class BitSetIndexHash {
	public static void main(String args[]){
		try{

			ArrayList<Fasta> genome=Fasta.readFasta(new File(args[0]));
			ArrayList<Fasta> sequences=Fasta.readFasta(new File(args[1]));
			File out=new File(args[2]);
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			ArrayList<BitSet> Genome=DNAmanipulations.FastaToBitSet(genome);
			int length=DNAmanipulations.codeDNA(sequences.get(0).getSequence()).length();
			BitSetIndexHash bsih=new BitSetIndexHash(Genome,length,true);
			System.out.println("Hash initialised...");
			for(int i=0;i<sequences.size();i++){
				System.out.println(sequences.get(i).ident);
				BitSet seq=DNAmanipulations.codeDNA(sequences.get(i).getSequence());
				if(length!=seq.length()){
					System.err.println(sequences.get(i).ident+" had not the correct length!");
					continue;
				}
				ArrayList<SequencePositions> al=bsih.getPos(seq);
				for(int j=0;j<al.size();j++){
					String orient=al.get(j).position>0?"f":"r";
					bw.write(sequences.get(i).ident+"\t"+genome.get(al.get(j).sequence).ident.replace("\t", "")+"\t"+orient+"\t"+Math.abs(al.get(j).position)+"\n");
				}
			}
			bw.close();
		}catch(IOException e){
			e.printStackTrace();
		}
	}
	
	HashMap<BitSet,ArrayList<SequencePositions>> index=new HashMap<BitSet, ArrayList<SequencePositions>>();
	HashMap<Integer,Integer> map;

	int sequence=0;
	boolean ori;
	public BitSetIndexHash(ArrayList<BitSet> input, int BitSetLength,boolean orientation){
		ori=orientation;

		for(int i=0;i<input.size();i++){
			indexing(input.get(i),BitSetLength);
		}
		
	}
	public BitSetIndexHash(ArrayList<BitSet> input, int BitSetLength,HashMap<Integer,Integer> pos,boolean orientation){
		map=pos;
		ori=orientation;

		for(int i=0;i<input.size();i++){
			indexing(input.get(i),BitSetLength);
		}
		
	}
	//length is supposed to denote the length of the bitset (bisetinstance.length())
	public BitSetIndexHash(BitSet input, int BitSetLength,boolean orientation){
		ori=orientation;
		indexing(input,BitSetLength);
	}
	
	

	private void indexing(BitSet input, int size){
		BitSet revInput=input;
		for(int i=0;i<input.length()-size;i+=2){
			BitSet key=input.get(i, i+size-1);
			key.set(size-1);
			SequencePositions temp=new SequencePositions(sequence,i/2);
			
			if(index.containsKey(key)){
				
				index.get(key).add(temp);
			}else{
				ArrayList<SequencePositions> al=new ArrayList<SequencePositions>();
				al.add(temp);
				index.put(key,al);
			}
			BitSet rev=revInput.get(i, i+size-1);
			rev.set(size-1);
			rev=DNAmanipulations.reverse(rev);
			int pos=ori?-i/2:i/2;
			temp=new SequencePositions(sequence,pos);
			if(index.containsKey(rev)){
				index.get(rev).add(temp);
			}else{
				ArrayList<SequencePositions> al=new ArrayList<SequencePositions>();
				al.add(temp);
				index.put(rev,al);
			}
		}
		sequence++;
	}
	public int getNumber(BitSet key){
		return index.get(key).size();
	}
	public ArrayList<SequencePositions> getPos(BitSet key){
		if(index.containsKey(key))
			return index.get(key);
		else return null;
	}
	//transforms sequence number and position into a continues sequence position (3 sequences a 30bp turns into one sequence of 90bp), where sequences on the 
	//opposite strain are kept negative
	public ArrayList<Integer> getPosMap(BitSet key){
		if(index.containsKey(key) && map!=null){
			ArrayList<SequencePositions> list=index.get(key);
			ArrayList<Integer> pos=new ArrayList<Integer>();
			for(int i=0;i<list.size();i++){
				int posi=list.get(i).position;
				
				int position=posi<0?-1*(map.get(list.get(i).sequence)+Math.abs(posi)):map.get(list.get(i).sequence)+Math.abs(posi);

				pos.add(position);
			}
			
			return pos;
		}
		else return null;
	}
	
	public int getNumber(ArrayList<BitSet> keys){
		int sum=0;
		for(int i=0;i<keys.size();i++){

			if(index.containsKey(keys.get(i)))
			sum+=index.get(keys.get(i)).size();
		}
		return sum;
	}
	
	public ArrayList<Integer> getPosMap(ArrayList<BitSet> keys){
		ArrayList<Integer> sum=new ArrayList<Integer>();
		for(int i=0;i<keys.size();i++){

			if(index.containsKey(keys.get(i))){
				ArrayList<SequencePositions> list=index.get(keys.get(i));
				ArrayList<Integer> pos=new ArrayList<Integer>();
				for(int j=0;j<list.size();j++){
					int posi=list.get(j).position;
					int position=posi<0?-1*(map.get(list.get(j).sequence)+Math.abs(posi)):map.get(list.get(j).sequence)+posi;
					pos.add(position);
				}
				sum.addAll(pos);
			}
		}
		return sum;
	}
	
	public ArrayList<SequencePositions> getPos(ArrayList<BitSet> keys){
		ArrayList<SequencePositions> sum=new ArrayList<SequencePositions>();
		for(int i=0;i<keys.size();i++){

			if(index.containsKey(keys.get(i))){
				
				sum.addAll(index.get(keys.get(i)));
			}
		}
		return sum;
	}
}
