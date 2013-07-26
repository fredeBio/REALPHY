package util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;

//import mini.ArtemisException;

public class ReadArtemis {
	//class to read artemis input files
	
	ArrayList<ArrayList<Feature>> artemisMap;
	
	public ReadArtemis(File artemis){
		artemisMap=new ArrayList<ArrayList<Feature>>();
		read(artemis);
	}
	
	public ArrayList<Feature> getFeatures(int pos1,int pos2){
		ArrayList<Feature> features=new ArrayList<Feature>();
		HashMap<Feature,Integer> check=new HashMap<Feature, Integer>();
		for(int i=pos1;i<=pos2;i++){
			for(int j=0;j<artemisMap.get(i).size();j++){
				if(!check.containsKey(artemisMap.get(i).get(j))){
					features.add(artemisMap.get(i).get(j));
					check.put(artemisMap.get(i).get(j),1);
				}
			}
		}
		return features;
	}
	
	public ArrayList<Feature> getFeatures(int pos1,int pos2,String featureName){
		ArrayList<Feature> features=new ArrayList<Feature>();
		HashMap<Feature,Integer> check=new HashMap<Feature, Integer>();
		for(int i=pos1;i<=pos2;i++){
			for(int j=0;j<artemisMap.get(i).size();j++){
				if(!check.containsKey(artemisMap.get(i).get(j)) && featureName.equals(artemisMap.get(i).get(j).name)){
					features.add(artemisMap.get(i).get(j));
					check.put(artemisMap.get(i).get(j),1);
				}
			}
		}
		return features;
	}
	
	public ArrayList<Feature> getFeatures(ArrayList<Integer> posList){
		ArrayList<Feature> features=new ArrayList<Feature>();
		HashMap<Feature,Integer> check=new HashMap<Feature, Integer>();
		for(int k=0;k<posList.size()-1;k+=2){
			int pos1 =posList.get(k);
			int pos2=posList.get(k+1);
			for(int i=pos1;i<=pos2;i++){
				for(int j=0;j<artemisMap.get(i).size();j++){
					if(!check.containsKey(artemisMap.get(i).get(j))){
						features.add(artemisMap.get(i).get(j));
						check.put(artemisMap.get(i).get(j),1);
					}
				}
			}
		}
		return features;
	}
	
	
	//true if feature is at this position in the genome
	public BitSet getBoolArray(String featureName){
		BitSet ba=new BitSet(artemisMap.size());
		ArrayList<Integer> posList=getPos(featureName);
		for(int j=0;j<posList.size()-1;j+=2){
			for(int k=posList.get(j);k<posList.get(j+1);k++){
				ba.set(k,true);
			}
		}
		ba.set(artemisMap.size()+1);	
		return ba;
		
	}
	public ArrayList<Integer> getPos(String featureName){
		ArrayList<Integer> posList=new ArrayList<Integer>();
		ArrayList<Feature> features=getFeatures(0, artemisMap.size()-1, featureName);
		for(int i=0;i<features.size();i++){
			ArrayList<Integer> start=features.get(i).start;
			ArrayList<Integer> end=features.get(i).end;
			for(int j=0;j<start.size();j++){
				posList.add(start.get(j));
				posList.add(end.get(j));
			}
		}
		return posList;
	}
	
	public static String featuresToString(ArrayList<Feature> features){
		StringBuilder feature=new StringBuilder("");
		for(int i=0;i<features.size();i++){
			Feature f=features.get(i);
			feature.append("FT   "+f.name+"            ");
			if(f.start.size()>1){
				feature.append("join(");
			}
			if(f.complement){
				feature.append("complement(");
			}
			for(int j=0;j<f.start.size()-1;j++){
				feature.append(f.start.get(j)+".."+f.end.get(j)+",");
			}
			feature.append(f.start.get(f.start.size()-1)+".."+f.end.get(f.end.size()-1));
			if(f.start.size()>1){
				feature.append(")");
			}
			if(f.complement){
				feature.append(")");
			}
			feature.append("\n");
			for(int j=0;j<f.properties.size();j++){
				feature.append("FT"+"   "+"                /"+f.properties.get(j).Name+"="+f.properties.get(j).Value+"\n");
			}
		}
		
		
		
		return feature.toString();
	}
	
	public void read(File artemis) {
		try{
			BufferedReader br = new BufferedReader(new FileReader(artemis));
			String line="";
			ArrayList<Property> properties=new ArrayList<Property>();
			int lineNumber=0;
			Feature f=new Feature();
			while((line=br.readLine())!=null){
				lineNumber++;
				if(lineNumber==1000){
					System.out.println("tst");
				}
				if(line.matches("FT   \\S+\\s+.+")){
					//feature & position
//					solve join & complement problem
					properties=new ArrayList<Property>(1);

					if(line.contains("complement(")){
						line=line.replace("complement(","");
						line=line.replace(")","");
						f.complement=true;
					}else{
						f.complement=false;
					}
					if(line.contains("join(")){
						line=line.replace("join(","");
						line=line.replace(")","");
					}
					if(line.contains("order(")){
						line=line.replace("order(","");
						line=line.replace(")","");
					}
					String[] features =line.split("\\.\\.>|\\s+<|\\s+|\\.\\.|,");
					f=new Feature();
					f.name=features[1];
					f.start=new ArrayList<Integer>(1);
					f.end=new ArrayList<Integer>(1);
					for(int i=2;i<features.length;i+=2){
						int pos1=Integer.parseInt(features[i]);
						int pos2=Integer.parseInt(features[i+1]);
						int t=pos1;
						pos1=pos1>pos2?pos2:pos1;
						pos2=pos1>pos2?t:pos2;

						if(pos2>=artemisMap.size()){
							for(int j=artemisMap.size();j<=pos2+1;j++){
								artemisMap.add(new ArrayList<Feature>(1));
							}
						}

						f.start.add(pos1);
						f.end.add(pos2);
						f.properties=properties;
						for(int j=pos1;j<=pos2;j++){
							artemisMap.get(j).add(f);
						}
					}
				}else if(line.matches("FT   \\s+/.+")){
					//new property
					String prop[]=line.split("\\s+|/|=");
					Property p=new Property();
					p.Name=prop[2];
					p.Value="";
					for(int i=3;i<prop.length-1;i++){
						p.Value+=prop[i]+" ";
					}
					p.Value+=prop[prop.length-1];
					properties.add(p);
				}else if(line.matches("FT   \\s+.+")){
					Property p=properties.get(properties.size()-1);
					String prop[]=line.split("\\s+");
					for(int i=3;i<prop.length-1;i++){
						p.Value+=prop[i]+" ";
					}

				}else{
				//	throw new ArtemisException("Error in artemis file \\"+artemis+" at line number "+lineNumber+".\n");
				}
				
			}

		}catch(IOException e){
			System.err.println(e.toString());
		}
	}
}

