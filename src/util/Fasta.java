package util;

import java.io.*;
import java.util.*;
import java.util.Map.Entry;


public class Fasta {
	String ident;
	private StringBuilder sequence;
	public static void main(String args[]){
		File fas=new File(args[0]);
		File outFolder=new File(args[1]);
		//int level=Integer.parseInt(args[2]);
		//writeIndividualFastaLeafs(Fasta.readFasta(fas), outFolder,level);
		//writePhylipLevel(Fasta.readFasta(fas),outFolder,level);
		writeIndividualFasta(Fasta.readFasta(fas), outFolder);
	}
	

	public static ArrayList<Fasta> sort(ArrayList<Fasta> fas,String[] order){
		HashMap<String,Fasta> fas2=Fasta.fasToFastaHash(fas,false);
		ArrayList<Fasta> newlist=new ArrayList<Fasta>();
		for(int i=0;i<order.length;i++){
			newlist.add(fas2.get(order[i]));
		}
		return newlist;
	}
	
	public static ArrayList<Fasta> sort(ArrayList<Fasta> order,ArrayList<Fasta> fas){
		HashMap<String,Fasta> fas2=Fasta.fasToFastaHash(fas,false);
		ArrayList<Fasta> newlist=new ArrayList<Fasta>();
		for(int i=0;i<order.size();i++){
			newlist.add(fas2.get(order.get(i).getIdent().trim()));
		}
		return newlist;
	}
	
	
	public static ArrayList<String> getSequences(ArrayList<Fasta> fas){
		ArrayList<String> seqs=new ArrayList<String>();
		for(int i=0;i<fas.size();i++){
			seqs.add(fas.get(i).getSequence());
		}
		return seqs;
	}
	
	public static ArrayList<String> getIdents(ArrayList<Fasta> fas){
		ArrayList<String> seqs=new ArrayList<String>();
		for(int i=0;i<fas.size();i++){
			seqs.add(fas.get(i).getIdent());
		}
		return seqs;
	}
	
	public static void writePhylipLevel(ArrayList<Fasta> fas,File out,int level){
		ArrayList<Fasta> temp=new ArrayList<Fasta>();
		for(int i=0;i<fas.size();i++){
			Fasta seq=fas.get(i);
			String ident=seq.getIdent().split("\\s+")[0];
			if(ident.length()==level+1){
				temp.add(seq);

			}
		}
		Fasta.writePhylip(temp, out,10);
	}

	
	public static void writeIndividualFasta(ArrayList<Fasta> fas,File outFolder){
		for(int i=0;i<fas.size();i++){
			Fasta seq=fas.get(i);
			ArrayList<Fasta> temp=new ArrayList<Fasta>();
			temp.add(seq);
			File out=new File(outFolder+"/"+seq.getIdent().split("\\s+")[0]+".fas");
			Fasta.write(temp, out);
		}
	}
	public static void writeIndividualFastaLeafs(ArrayList<Fasta> fas,File outFolder,int level){
		for(int i=0;i<fas.size();i++){
			Fasta seq=fas.get(i);
			String ident=seq.getIdent().split("\\s+")[0];
			if(ident.length()==level+1){
				ArrayList<Fasta> temp=new ArrayList<Fasta>();
				seq.setIdent(ident);
				temp.add(seq);

				File out=new File(outFolder+"/"+ident+".fas");
				Fasta.write(temp, out);
			}
		}
	}
	
	public static void write(HashMap<String,String> fas,File out){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			Iterator<Entry<String,String>> it=fas.entrySet().iterator();
			while(it.hasNext()){
				Entry<String,String> e=it.next();
				bw.write(">"+e.getKey()+"\n"+e.getValue()+"\n");
			}
			
			bw.close();
		}catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	public void setIdent(String id){
		ident=id;
	}
	public String getSequence(){
		return sequence.toString();
	}
	public Fasta(String Ident,String Sequence){
		ident=Ident;
		sequence=new StringBuilder(Sequence);
	}
	public String toString(){
		String dummy=">"+ident+"\n"+sequence+"\n";
		return sequence.length()>0?dummy:"";
	}
	public static Fasta makeFasta(Info interval,String genome,boolean translate){
		int start=interval.getStart();
		int end=interval.getEnd();
		String seq=genome.substring(start-1, end);
		String name=interval.info+"_"+start+"_"+end;
		if(interval.info.endsWith("complement")){
			seq=DNAmanipulations.reverse(seq);
		}
		if(translate){
			seq=DNAmanipulations.translate(seq,DNAmanipulations.code());
		}
		 
	return new Fasta(name,seq);
	}
	
	public void setSequence(String Sequence){
		sequence=new StringBuilder(Sequence);
	}
	
	public static void write(ArrayList<Fasta> fasta,File out){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			for(int i=0;i<fasta.size();i++){
				bw.write(fasta.get(i)+"");
			}
			bw.close();
		}catch(IOException e){
			e.printStackTrace();
		}
	}
	public static void write(ArrayList<Fasta> fasta,File out,boolean append){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(out,append));
			for(int i=0;i<fasta.size();i++){
				bw.write(fasta.get(i)+"");
			}
			bw.close();
		}catch(IOException e){
			e.printStackTrace();
		}
	}
	public static String getSequence(Info interval,String genome,boolean translate){
			int start=interval.getStart();
			int end=interval.getEnd();
			if(start>end){
				System.err.println("Warning, start position ("+start+") greater than end position ("+end+").");
				return "";
			}
			String seq=genome.substring(start-1, end);
			if(interval.info.endsWith("complement")){
				seq=DNAmanipulations.reverse(seq);
			}
			if(translate){
				seq=DNAmanipulations.translate(seq,DNAmanipulations.code());
			}
			 
		return seq;
	}
	
	public static String makeIdent(int idLength,String fasIdent){
		String ident=fasIdent.split("\\s+")[0];
		if(ident.length()>idLength-1){
			ident=ident.substring(0,idLength-1);
			ident+=" ";
		}else{
			int length=ident.length();
			for(int j=0;j<idLength-length;j++){
				ident+=" ";
			}
		}
		return ident;
	}
	
	public static String makeIdentNumbered(int idLength,String fasIdent,int number){
		String ident=fasIdent.split("\\s+")[0];
		int digits=number==0?1:(int)Math.log10(number)+1;
		if(ident.length()>idLength-1-digits){
			ident=ident.substring(0,idLength-1-digits);
			ident+=number+" ";
		}else{
			int length=ident.length()+digits;
			ident+=number;
			for(int j=0;j<idLength-length;j++){
				ident+=" ";
			}
		}
		return ident;
	}
	
	public static int writePhylip(ArrayList<Fasta> fas,File out,int idLength){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			int seqLength=0;
			for(int i=0;i<fas.size();i++){
				if(i==0){
					seqLength=fas.get(i).getSequence().length();
					bw.write("\t"+fas.size()+" "+seqLength+"\n");
				}
				String ident=makeIdent(idLength,fas.get(i).getIdent());
				bw.write(ident+fas.get(i).getSequence()+"\n");
			}
			
			bw.close();
			return seqLength;
		}catch(IOException e){
			e.printStackTrace();
		}
		return -1;
	}

	public static int writePhylipNumbered(ArrayList<Fasta> fas,File out,int idLength){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			int seqLength=0;
			int size=fas.size();
			for(int i=0;i<size;i++){
				if(i==0){
					seqLength=fas.get(i).getSequence().length();
					bw.write("\t"+fas.size()+" "+seqLength+"\n");
				}
				String ident=makeIdentNumbered(idLength,fas.get(i).getIdent(),i);
				bw.write(ident+fas.get(i).getSequence()+"\n");
			}
			
			bw.close();
			return seqLength;
		}catch(IOException e){
			e.printStackTrace();
		}
		return -1;
	}

	
	public static ArrayList<Fasta> readPhylipRAxML(File in){
		ArrayList<Fasta> fas=new ArrayList<Fasta>();
		HashMap<String,StringBuffer> hm=readPhylipHM(in);
		fas=HashToFas(hm,hm.keySet().toArray(new String[0]));
		
		return fas;
	}
	
	public static ArrayList<Fasta> readPhylip(File in,int idLength){
		ArrayList<Fasta> fas=new ArrayList<Fasta>();
		HashMap<String,StringBuffer> hm=readPhylipHM(in,idLength);
		fas=HashToFas(hm,hm.keySet().toArray(new String[0]));
		System.err.println("Problem? readphylip idlength "+in);
		return fas;
	}

	public static HashMap<String,StringBuffer> readPhylipHM(File in,int idLength){
		HashMap<String,StringBuffer> hm=new HashMap<String, StringBuffer>();

		try{
			BufferedReader br=new BufferedReader(new FileReader(in));
			String line="";
			int i=0;
			ArrayList<String> idents=new ArrayList<String>();
			while((line=br.readLine())!=null){
				if(i==0){
					i++;
				}else{
					if(!line.matches("\\s+")){
						String id=line.substring(0,idLength);
						String seq=line.substring(idLength).trim().toUpperCase();
						//System.out.println(id);
						if(hm.containsKey(id)){
							hm.get(id).append(seq);
						}else{
							idents.add(id);
							hm.put(id,new StringBuffer(seq));
						}
					}
				}
				
			}
			br.close();	
		}catch(IOException e){
			e.printStackTrace();
		}
		return hm;
	}
	
	public static HashMap<String,StringBuffer> readPhylipHM(File in){
		HashMap<String,StringBuffer> hm=new HashMap<String, StringBuffer>();

		try{
			BufferedReader br=new BufferedReader(new FileReader(in));
			String line="";
			int i=0;
			ArrayList<String> idents=new ArrayList<String>();
			while((line=br.readLine())!=null){
				if(i==0){
					i++;
				}else{
					if(!line.matches("\\s+")){
						String split[]=line.split("\\s+");
						String id=split[0];
						
						//System.out.println(id);
						if(hm.containsKey(id)){
							hm.get(id).append(split[split.length-1].trim().toUpperCase());
						}else{
							idents.add(id);
							hm.put(id,new StringBuffer(split[split.length-1].trim().toUpperCase()));
						}
					}
				}
				
			}
			br.close();	
		}catch(IOException e){
			e.printStackTrace();
		}
		return hm;
	}
	
	public static ArrayList<Fasta> readPhylip(File in,String[] order){
		ArrayList<Fasta> fas=new ArrayList<Fasta>();
			HashMap<String,StringBuffer> hm=readPhylipHM(in);
			fas=HashToFas(hm,order);
			return fas;
	}
	
	public static ArrayList<Fasta> readPhylip(File in){
		ArrayList<Fasta> fas=new ArrayList<Fasta>();
		HashMap<String,StringBuffer> hm=readPhylipHM(in);
		fas=HashToFas(hm,hm.keySet().toArray(new String[0]));
		System.err.println("Problem readPhylip? "+in);
		return fas;
	}
	
	public static String getColumn(ArrayList<Fasta> fas,int col){
		StringBuffer Col=new StringBuffer();
		for(int i=0;i<fas.size();i++){
			Col.append(fas.get(i).getSequence().charAt(col));
		}
		return Col.toString();
	}
	
	public static ArrayList<Fasta> HashToFas(HashMap<String,StringBuffer> hm,String[] idents){
		ArrayList<Fasta> fas=new ArrayList<Fasta>();
		for(int i=0;i<idents.length;i++){
			
			fas.add(new Fasta(idents[i],hm.get(idents[i]).toString()));
			
		}
		return fas;
	}
	
	public static HashMap<String,String> fasToHash(ArrayList<Fasta> fas,boolean resolveGI){
		HashMap<String,String> hash=new HashMap<String, String>();
		for(int i=0;i<fas.size();i++){
			String ident=fas.get(i).getIdent();
			if(!resolveGI||!ident.startsWith("gi|")){
				String[] split=ident.split("\\s+");
				hash.put(split[0], fas.get(i).getSequence());
			}else{
				String[] split=ident.split("\\||\\.");
				hash.put(split[3], fas.get(i).getSequence());
			}
		}
		return hash;
	}
	
	public static HashMap<String,Fasta> fasToFastaHash(ArrayList<Fasta> fas,boolean resolveGI){
		HashMap<String,Fasta> hash=new HashMap<String, Fasta>();
		for(int i=0;i<fas.size();i++){
			String ident=fas.get(i).getIdent();
			if(!resolveGI||!ident.startsWith("gi|")){
				String[] split=ident.split("\\s+");
				hash.put(split[0], fas.get(i));
			}else{
				String[] split=ident.split("\\||\\.");
				hash.put(split[3], fas.get(i));
			}
		}
		return hash;
	}
	
	
	public static ArrayList<Fasta> readFasta(File f){
		try{
			BufferedReader sequence=new BufferedReader(new FileReader(f));
			return read(sequence);
		}catch(IOException e){
			System.err.println(e.toString());
		}
		return null;
	}
	
	public static ArrayList<Fasta> readFasta(File f,HashMap<String,Integer> order){
		try{
			
			BufferedReader sequence=new BufferedReader(new FileReader(f));
			return read(sequence,order);
		}catch(IOException e){
			System.err.println(e.toString());
		}
		return null;
	}
	
	private static ArrayList<Fasta> read(BufferedReader br){
		ArrayList<Fasta> fasta=new ArrayList<Fasta>();
		try{
			String line="";
			while((line=br.readLine())!=null){
				if(line.startsWith(">")){
					fasta.add(new Fasta(line.replace(">",""),""));
				}else{
					line=line.replaceAll("\\s", "");
					if(fasta.size()>0){
						fasta.get(fasta.size()-1).sequence.append(line.toUpperCase());
					}
				}
			}

		}catch(IOException e){
			System.err.println(e.toString());
		}
		return fasta;
	}
	
	private static ArrayList<Fasta> init(int size){
		ArrayList<Fasta> al=new ArrayList<Fasta>();
		for(int i=0;i<size;i++){
			al.add(new Fasta("",""));
		}
		return al;
	}
	

	
	private static ArrayList<Fasta> read(BufferedReader br,HashMap<String,Integer> order){
		ArrayList<Fasta> fasta=init(order.size());
		try{
			String line="";
			int pos=0;
			while((line=br.readLine())!=null){
				if(line.startsWith(">")){
					String ident=line.replace(">","");
					if(order.containsKey(ident)){
						pos=order.get(ident);
						fasta.set(pos,new Fasta(ident,""));
					}else{
						System.err.println(ident+" was not found in order file!");
						System.exit(-1);
					}
				}else{
					line.replaceAll("\\s", "");
					if(fasta.size()>0){
						fasta.get(pos).sequence.append(line.toUpperCase());
					}
				}
			}

		}catch(IOException e){
			System.err.println(e.toString());
		}
		return fasta;
	}

	
	public String getIdent(){
		return ident;
	}
	public static ArrayList<Fasta> translate(ArrayList<Fasta> seq){
		ArrayList<Fasta> trans=new ArrayList<Fasta>();
		for(int i=0;i<seq.size();i++){
			trans.add(new Fasta(seq.get(i).ident,DNAmanipulations.translate(seq.get(i).getSequence(), DNAmanipulations.code())));
		}
		return trans;
	}
	
	public static ArrayList<Fasta> makeFasta(ArrayList<Info> intervals,String genome,boolean translate){
		genome=genome.toUpperCase();
		ArrayList<Fasta> list=new ArrayList<Fasta>();
		for(int i=0;i<intervals.size();i++){
			int start=intervals.get(i).getStart();
			int end=intervals.get(i).getEnd();
			String seq=genome.substring(start-1, end);
			String name=intervals.get(i).info+"_"+start+"_"+end;
			if(intervals.get(i).info.endsWith("complement")){
				seq=DNAmanipulations.reverse(seq);
			}
			if(translate){
				seq=DNAmanipulations.translate(seq,DNAmanipulations.code());
			}
			list.add(new Fasta(name,seq));
			 
		}
		return list;
	}
	
}
