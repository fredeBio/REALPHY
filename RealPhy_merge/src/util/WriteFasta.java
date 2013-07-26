package util;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map.Entry;


public class WriteFasta {
	public static void write(HashMap<String,String> fasta,File out){
		try{
			Iterator<Entry<String,String>> it=fasta.entrySet().iterator();
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			while(it.hasNext()){
				Entry<String,String> e=it.next();
				bw.write(">"+e.getKey()+"\n"+e.getValue()+"\n");
			}
			bw.close();
		}catch(IOException e){
			System.err.println(e.toString());
		}
	}


}
