package util;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

public class WriteArtemis {
	public static void write(ArrayList<Info> info,File art){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(art));
			for(int i=0;i<info.size();i++){
				String ide=info.get(i).info.length()>5?info.get(i).info.substring(0,5):info.get(i).info;
				StringBuilder sb=new StringBuilder(ide);
				for(int j=0;j<5-ide.length();j++){
					sb.append(" ");
				}
				ide=sb.toString();
				int start=info.get(i).getStart();
				int end=info.get(i).getEnd();
				String inf=info.get(i).info;
				if(start<end)
					bw.write("FT   "+ide+"             "+start+".."+end+"\n");
				else bw.write("FT   "+ide+"             "+"complement("+start+".."+end+")"+"\n");
				bw.write("FT                   /note=\""+inf+"\"\n");
			}
			bw.close();
		}catch(IOException e){
			System.err.println(e.toString());
		}
	}
	public static void write(ArrayList<Integer> start,ArrayList<Integer> end,ArrayList<String> id,File art){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(art));
			for(int i=0;i<start.size();i++){
				String ide=id.get(i).length()>5?id.get(i).substring(0,5):id.get(i);
				StringBuilder sb=new StringBuilder(ide);
				for(int j=0;j<5-ide.length();j++){
					sb.append(" ");
				}
				ide=sb.toString();
				if(start.get(i)<end.get(i))
					bw.write("FT   "+ide+"             "+start.get(i)+".."+end.get(i)+"\n");
				else bw.write("FT   "+ide+"             "+"complement("+start.get(i)+".."+end.get(i)+")"+"\n");
			}
			bw.close();
		}catch(IOException e){
			System.err.println(e.toString());
		}
	}

}
