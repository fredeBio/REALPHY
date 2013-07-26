package util;

import java.io.*;
import java.util.*;
import java.util.regex.*;




public class ReadGenbank {
	
	ArrayList<Fasta> sequence=new ArrayList<Fasta>();
	
	HashMap<String,InfoTree> infoHash=new HashMap<String, InfoTree>();
	ArrayList<String> ids=new ArrayList<String>();
	public ArrayList<Info> getFeatures(Info pos,String id){
		ArrayList<Info> al=new ArrayList<Info>();
		infoHash.get(id).search(pos,al);
	//	System.out.println(infoTree.find(new Info(1090939,1092767,"test")).info);
		return al;
	}
	public ArrayList<Info> getFeatures(Info pos){
		ArrayList<Info> al=new ArrayList<Info>();
		infoHash.get(ids.get(0)).search(pos,al);
	//	System.out.println(infoTree.find(new Info(1090939,1092767,"test")).info);
		return al;
	}
	
	public  ArrayList<String> getIds(){
		return ids;
	}
	public static void main(String args[]){
		ReadGenbank rgb=new ReadGenbank(new File(args[0]));
		System.out.println(rgb.infoHash.get("pfs"));
		
	}
	
	public ReadGenbank(File GBK,String... features){
		parseGBK(GBK,features);
	}
	
	public InfoTree getInfoTree(String feature,String id){
		ArrayList<Info> features=getIntervals(feature,id);
		InfoTree itree=new InfoTree();
		for(int i=0;i<features.size();i++){
			itree.insert(features.get(i));
		}
		return itree;
	}
	public InfoTree getInfoTree(String feature){
		ArrayList<Info> features=getIntervals(feature,ids.get(0));
		InfoTree itree=new InfoTree();
		for(int i=0;i<features.size();i++){
			itree.insert(features.get(i));
		}
		return itree;
	}
	public InfoTree getWholeInfoTree(){
		return infoHash.get(ids.get(0));
	}
	public InfoTree getWholeInfoTree(String id){
		return infoHash.get(id);
	}
	
	public ArrayList<Info> getIntervals(String feature,String id){
		ArrayList<Info> al=new ArrayList<Info>();
		al=infoHash.get(id).parseTree(feature);
		return al;
	}
	public ArrayList<Info> getIntervals(HashMap<String,Boolean> feature){
		ArrayList<Info> al=new ArrayList<Info>();
		al=infoHash.get(ids.get(0)).parseTree(feature);
		return al;
	}
	public ArrayList<Info> getIntervals(String feature){
		ArrayList<Info> al=new ArrayList<Info>();
		al=infoHash.get(ids.get(0)).parseTree(feature);
		return al;
	}
	public ArrayList<Info> getIntervals(String feature,int minLength,String id){
		ArrayList<Info> al=new ArrayList<Info>();
		al=infoHash.get(id).parseTree(feature,minLength);
		return al;
	}
	public ArrayList<Info> getIntervals(String feature,int minLength){
		ArrayList<Info> al=new ArrayList<Info>();
		al=infoHash.get(ids.get(0)).parseTree(feature,minLength);
		return al;
	}
	public ArrayList<Info> getFeatures(String id){
		return infoHash.get(id).parseTree();
	}
	
	public ArrayList<Fasta> getFeaturesSequences(){
		ArrayList<Fasta> features=new ArrayList<Fasta>();
		int count=0;
		for(int i=0;i<sequence.size();i++){
			ArrayList<Info> infList=infoHash.get(sequence.get(i).getIdent()).parseTree();
			String wholeSeq=sequence.get(i).getSequence();
			for(int j=0;j<infList.size();j++){
				Info inf=infList.get(j);
				
				String ident=count+" "+inf.info;
				String sequence=inf.orient=='+'?wholeSeq.substring(inf.start, inf.end):DNAmanipulations.reverse(wholeSeq.substring(inf.start, inf.end));
				features.add(new Fasta(ident,sequence));
				count++;
			}
		}
		return features;
	}
	
	public ArrayList<Info> getFirstSequenceFeatures(){
		return infoHash.get(ids.get(0)).parseTree();
	}
	
	public HashMap<String,ArrayList<Info>> getFeatures(){
		HashMap<String,ArrayList<Info>> hm=new HashMap<String, ArrayList<Info>>();
		for(int i=0;i<ids.size();i++){
			hm.put(ids.get(i), infoHash.get(ids.get(i)).parseTree());
		}
		return hm;
	}
	
	public ArrayList<Fasta> getSequence(){
		return sequence;
	}
	private static HashMap<String,Boolean> makeHash(String... args){
		HashMap<String,Boolean> hm=new HashMap<String, Boolean>();
		for(int i=0;i<args.length;i++){
			hm.put(args[i],true);
		}
		return hm;
	}
	private void parseGBK(File GBK,String... onlyFeature){
		//System.out.println(GBK);
		try{
			HashMap<String,Boolean> oFeature=onlyFeature.length>0?makeHash(onlyFeature):null;
			BufferedReader br=new BufferedReader(new FileReader(GBK));
			String line="";
			int start=-1;
			int end=-1;
			String info="";
			boolean found=false;
			boolean complement=false;
			String locus="";
			InfoTree infoTree=new InfoTree();
			boolean origin=false;
			StringBuffer tempSeq=new StringBuffer();
			int codon_start=-1;
			boolean pseudo=false;
			String feature="";
			while((line=br.readLine())!=null){
				if(line.matches("^\\s+\\S+\\s+complement\\(join.+")||line.matches("^\\s+\\S+\\s+join\\(.+")){
					found=false;
					if(info.length()>0 && start>-1){
						info=complement?info+" complement":info;
						char orient=complement?'-':'+';
						if(oFeature==null||oFeature.containsKey(feature))infoTree.insert((new Info(start,end,info,orient)).setPseudo(pseudo).setFeature(feature));
						info=new String("");
						pseudo = false;
					}
				}else if(line.startsWith("LOCUS")){
					if(start>-1){
						if(info.length()>0){
							info=complement?info+" complement":info;
							char orient=complement?'-':'+';
							if(oFeature==null||oFeature.containsKey(feature))infoTree.insert((new Info(start,end,info,orient)).setPseudo(pseudo).setFeature(feature));
							infoHash.put(locus, infoTree);
							ids.add(locus);
							pseudo=false;
						}		
					}
					origin=false;
					locus=line.split("\\s+")[1];
					infoTree=new InfoTree();
					start=-1;
				}else if((line.matches("^\\s+\\S+\\s+complement\\(\\d+.*")||line.matches("^\\s+\\S+\\s+<*\\d+\\.\\.>*\\d+\\s*") )){
					Matcher m=Pattern.compile("^\\s+(\\S+)\\s+complement\\(<*(\\d+)\\.\\.>*(\\d+)\\)").matcher(line);
					if(info.length()>0 && start>-1){
						info=complement?info+" complement":info;
						char orient=complement?'-':'+';
						if(oFeature==null||oFeature.containsKey(feature)){
							if(codon_start>-1){
								Info inf=new Info(start,end,info).setOrient(orient).setCodonStart(codon_start).setPseudo(pseudo).setFeature(feature);
								infoTree.insert(inf);
								System.out.println(feature+ " sss "+inf.getFeature());
							}
							else {
								infoTree.insert(new Info(start,end,info).setOrient(orient).setPseudo(pseudo).setFeature(feature));
							}
						}
						pseudo=false;
						//System.out.println(info);
					}
					complement=true;
					if(!m.find()){
						m=Pattern.compile("^\\s+(\\S+)\\s+<*(\\d+)\\.\\.>*(\\d+)\\s*").matcher(line);
						complement=false;
						if(!m.find()){
							m=Pattern.compile("^\\s+(\\S+)\\s+<*>*(\\d+)").matcher(line);
							complement=false;
							if(!m.find()){
								m=Pattern.compile("^\\s+(\\S+)\\s+complement\\(<*>*(\\d+)\\)").matcher(line);
								complement=true;
							}
						}
					}
					m.reset();
					m.find();
					start=Integer.parseInt(m.group(2));
					
					if(m.groupCount()==3)end=Integer.parseInt(m.group(3));
					else end=start;
					info=" ";
					feature=m.group(1);
					found=true;
				}else if(found && line.matches("\\s+/note=\".+")){
					Matcher m=Pattern.compile("\\s+/(note=\".+)").matcher(line);
					m.find();
					info+=m.group(1)+" ";
					info=info.replace('"', ' ');

				}else if(found && line.matches("\\s+/product=\".+")){
					Matcher m=Pattern.compile("\\s+/(product=\".+)").matcher(line);
					m.find();
					info+=m.group(1).replace(' ', '_')+" ";
					info=info.replace('"', ' ');

				}else if(found && line.matches("\\s+/colour=.+")){
					Matcher m=Pattern.compile("\\s+/colour=(.+)").matcher(line);
					m.find();
					info+="colour= "+m.group(1)+" ";
					info=info.replace('"', ' ');

				}else if(found && line.matches("\\s+/EC_number=\".+")){
					Matcher m=Pattern.compile("\\s+/(EC_number=\".+)").matcher(line);
					m.find(); 
					info+=m.group(1)+" ";
					info=info.replace('"', ' ');

				}else if(found && line.matches("\\s+/codon_start=\".+")){
					Matcher m=Pattern.compile("\\s+/codon_start=\"(.+)").matcher(line);
					m.find();
					codon_start=Integer.parseInt(m.group(1));
				}else if(found && line.matches("\\s+/locus_tag=\".+")){
					Matcher m=Pattern.compile("\\s+/(locus_tag=\".+)").matcher(line);
					m.find();
					info+=m.group(1)+" ";
					info=info.replace('"', ' ');

				}else  if(found && line.matches("\\s+/gene=\".+") ){
					Matcher m=Pattern.compile("\\s+/(gene=\".+)").matcher(line);
					m.find();
					info+=m.group(1)+" ";
					info=info.replace('"', ' ');

				}else if(found && line.matches("\\s+/mobile_element=\".+") ){
					Matcher m=Pattern.compile("\\s+/(mobile_element=\".+)").matcher(line);
					m.find();
					info+=m.group(1)+" ";
					info=info.replace('"', ' ');

				}else if(found && line.matches("\\s+/pseudo\\s+") ){
					pseudo=true;
				}else if(line.startsWith("//")){
					origin=false;
					sequence.add(new Fasta(locus,tempSeq.toString()));
					tempSeq=new StringBuffer();
				}else if(line.startsWith("ORIGIN")){
					origin=true;
				}else if(origin){
					String split[]=line.split("\\s+");
					for(int i=2;i<split.length;i++)
						tempSeq.append(split[i]);
				}
			}
			if(found){
				info=complement?info+" complement":info;
				char orient=complement?'-':'+';
				if(oFeature==null||oFeature.containsKey(feature)){
					if(codon_start>-1)infoTree.insert((new Info(start,end,info)).setOrient(orient).setCodonStart(codon_start).setPseudo(pseudo).setFeature(feature));
					else infoTree.insert((new Info(start,end,info)).setOrient(orient).setPseudo(pseudo).setFeature(feature));
				}
				
				infoHash.put(locus,infoTree);

				ids.add(locus);
			}
			br.close();
		}catch(IOException e){
			e.printStackTrace();
		}		
	}






	
}
