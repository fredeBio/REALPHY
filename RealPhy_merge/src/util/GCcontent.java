package util;


import java.io.File;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map.Entry;

//measures the GC content of the sequences in a fasta file and the overall GC content
//input fasta file

public class GCcontent {
	public static void main(String args[]){
		File fasta = new File(args[0]);
		HashMap<String,StringBuilder> fas=ReadFasta.readFasta(fasta);
		Iterator<Entry<String,StringBuilder>> it=fas.entrySet().iterator();
		HashMap<Character,Integer> overall=new HashMap<Character, Integer>();
		while(it.hasNext()){
			Entry<String,StringBuilder> e =it.next();
			HashMap<Character,Integer> current=countATGC(e.getValue().toString());
			output(current);
			System.out.println("Sequencelength: "+e.getValue().length());
			overall=add(overall,current);
		}
		System.out.println("Overall:");
		output(overall);
	}
	public static HashMap<Character,Integer> countATGC(String sequence){
		HashMap<Character,Integer> atgc=new HashMap<Character, Integer>();
		atgc.put('A', 0);
		atgc.put('T', 0);
		atgc.put('G', 0);
		atgc.put('C', 0);

		for(int i = 0;i<sequence.length()/2;i++){
			Character c=Character.toUpperCase(sequence.charAt(i));
			if(atgc.containsKey(c))atgc.put(c, atgc.get(c)+1);
			
		}
		
		return atgc;
	}
	
	public static HashMap<Character,Integer> add(HashMap<Character,Integer> h1,HashMap<Character,Integer> h2){
		Iterator<Entry<Character,Integer>> it=h1.entrySet().iterator();
		while(it.hasNext()){
			Entry <Character,Integer> e=it.next();
			if(h2.containsKey(e.getKey())){
				h2.put(e.getKey(),e.getValue()+h2.get(e.getKey()));
			}else{
				h2.put(e.getKey(), e.getValue());
			}
		}
		
		return h2;
	}
	
	public static void output(HashMap<Character,Integer> h){
		double gcCont=0;
		gcCont=(h.get('G')+h.get('C')+0.0)/(h.get('C')+h.get('G')+h.get('A')+h.get('T')+0.0)*100;
		
		System.out.println("G: "+h.get('G')+"\nC: "+h.get('C')+"\nA: "+h.get('A')+"\nT: "+h.get('T')+"\nGC content:"+gcCont+"%");
	}
}
