package util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;

public class ReadFasta {
	
	private static HashMap<String, BitSet> readBitSet(BufferedReader br){
		HashMap<String, BitSet> fasta=new HashMap<String, BitSet>();
		String line="";
		String id="";
		try{
			int size=0;
			int k=1;
		while((line=br.readLine())!=null){
			if(size>k*1000000  ){
				
				System.out.println(k +"mio bases read...");
				k++;
			}
			if(line.startsWith(">")){
				id=line;
			}else
				if(id.length()==0){
					continue;
				}else{
					line.replaceAll("\\s", "");
					size+=line.length();
					if(line.contains("N")){
						System.err.println("Ns were converted to Cs!");
						line=line.toUpperCase().replace('N','C');
					}
					if(fasta.containsKey(id)){
						fasta.put(id,DNAmanipulations.append(fasta.get(id),DNAmanipulations.codeDNA(line.toUpperCase())));
					}else{
						fasta.put(id,DNAmanipulations.codeDNA(line.toUpperCase()));
					}
				}
			
		}
		}catch(IOException e){
			e.printStackTrace();
		}
		return fasta;
	}
	
	public static HashMap<String, BitSet> readFastaBitSet(File f){
		try{
			BufferedReader br=new BufferedReader(new FileReader(f));
			return readBitSet(br);
		}catch(IOException e){
			System.err.println(e.toString());
		}
		return null;
	}
	public static HashMap<String, StringBuilder> readFasta(File f){
		try{
			BufferedReader sequence=new BufferedReader(new FileReader(f));
			return read(sequence);
		}catch(IOException e){
			System.err.println(e.toString());
		}
		return null;
	}

	private static HashMap<String,StringBuilder> read(BufferedReader br){
		HashMap<String, StringBuilder> fasta=new HashMap<String, StringBuilder>();
		try{
			String line="";
			String id="";
			while((line=br.readLine())!=null){
				if(line.startsWith(">")){
					//changed from id=line;
					if(!line.startsWith(">gi|"))id=line.split("\\s+")[0];
					else id=line.split("\\s+|\\||\\.")[3];
				}else
					if(id.length()==0){
						continue;
					}else{
						line.replaceAll("\\s", "");
						if(fasta.containsKey(id)){
							fasta.get(id).append(line);
						}else{
							fasta.put(id,new StringBuilder(line));
						}
					}
				
			}
		}catch(IOException e){
			System.err.println(e.toString());
		}
		return fasta;
	}
	
	public static HashMap<String, StringBuilder> readFasta(String seq){

			BufferedReader sequence=new BufferedReader(new StringReader(seq));


		return read(sequence);
	}
	
	public static ArrayList<Integer> getPositions(String Sequence,String Genome){
		ArrayList<Integer> al=new ArrayList<Integer>();
		int New=0;
		int Old=0;
		Sequence=Sequence.toUpperCase();
		Genome=Genome.toUpperCase();
		while((New=Genome.indexOf(Sequence,Old))!=-1){
			al.add(New);
			Old=New+1;
		}
		
		return al;
	}
	public static ArrayList<Integer> getPositionsBoth(String Sequence,String Genome){
		ArrayList<Integer> al=new ArrayList<Integer>();
		int Old=0;
		Sequence=Sequence.toUpperCase();
		Genome=Genome.toUpperCase();
		String ReverseSequence=DNAmanipulations.reverse(Sequence);
		boolean end=Genome.indexOf(Sequence,Old)==-1 && Genome.indexOf(ReverseSequence,Old)==-1?true:false;
		while(!end){
			int seqfind=Genome.indexOf(Sequence,Old);
			int seqfindRev=Genome.indexOf(ReverseSequence,Old);
			if(seqfind >-1 && seqfindRev >-1){
				if(seqfind<seqfindRev){
					al.add(seqfind+1);al.add(seqfind+Sequence.length());

					Old=seqfind+1;
				}else{
					al.add(seqfindRev+Sequence.length());al.add(seqfindRev+1);

					Old=seqfindRev+1;
				}
			}else{
				if (seqfind==-1 &&seqfindRev==-1){
					end=true;

					return al;
				}else{
					if(seqfind>-1){
						al.add(seqfind+1);al.add(seqfind+Sequence.length());
						Old=seqfind+1;
					}else{
						al.add(seqfindRev+Sequence.length());al.add(seqfindRev+1);

						Old=seqfindRev+1;
					}
				}
				
			}
			
		}

		return al;
	}
	
	BufferedReader br;
	String previousID=new String();
	String previousLine=new String();
	int bases;
	int w;
	
	//TODO: Just works when line length is longer than search length...and when line length is smaller than bases...
	
	public HashMap<String,BitSet> read(){
		int i=0;
		String line;
		String id=previousID;
		HashMap<String,BitSet> fasta=new HashMap<String, BitSet>();
		boolean start=true;
		try{
			while(i<bases && (line=br.readLine())!=null){
				if(line.startsWith(">")){
					id=line;
					previousID=id;
					if(start){
						previousLine=new String();
					}
				}else

					if(id.length()==0){
						continue;
					}else{
						
						line.replaceAll("\\s", "");
						if(line.toUpperCase().contains("N")){
							line=line.toUpperCase().replace('N', 'C');
							System.err.println("WARNING: Ns replaced by Cs!");
						}
						if(start && previousLine.length()>0){
							String p=previousLine.substring(previousLine.length()-w+1);
							i+=p.length();
							if(fasta.containsKey(id)){
								fasta.put(id,DNAmanipulations.append(fasta.get(id),DNAmanipulations.codeDNA(p.toUpperCase())));
							}else{
								fasta.put(id,DNAmanipulations.codeDNA(p.toUpperCase()));
							}
							
						}
						start=false;
						i+=line.length();
						previousLine=line;
						if(fasta.containsKey(id)){
							fasta.put(id,DNAmanipulations.append(fasta.get(id),DNAmanipulations.codeDNA(line.toUpperCase())));
						}else{
							fasta.put(id,DNAmanipulations.codeDNA(line.toUpperCase()));
						}
					}
			}
		}catch(IOException e){
			System.err.println(e.toString());
		}
		return fasta;

	}
	
	public ReadFasta(File genome,int basepairs,int wordsize){
		bases=basepairs;
		w=wordsize;
		try{
			br=new BufferedReader(new FileReader(genome));
		}catch(IOException e){
			System.err.println(e.toString());
		}
	}
	public ReadFasta(String genome,int basepairs,int wordsize){
		bases=basepairs;
		w=wordsize;
			br=new BufferedReader(new StringReader(genome));

	}
	
	public static HashMap<String, BitSet> readFastaBitSet(String sequence) {
		BufferedReader br=new BufferedReader(new StringReader(sequence));
		return readBitSet(br);
	}
}
