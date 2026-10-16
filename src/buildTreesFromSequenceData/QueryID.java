package buildTreesFromSequenceData;

import java.io.Serializable;
import java.util.*;

public class QueryID implements Serializable{
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	
	//Hash<contig number,Hash<position,ID>
	HashMap<Integer,HashMap<Integer,Integer>> qID=new HashMap<Integer, HashMap<Integer,Integer>>();
	ArrayList<Query> qIDList=new ArrayList<Query>();
	
	ArrayList<String> contigsList=new ArrayList<String>();
	HashMap<String,Integer> contigs=new HashMap<String, Integer>();
	
	HashMap<Integer,Boolean> removed=new HashMap<Integer, Boolean>();

	
	public void remove(int key){
		removed.put(key,true);
	}
	
	
	
	public int get(String query){
		Query q=parse(query);
		int pos=q.pos;
		int contig=q.contig;
		if(qID.containsKey(contig)){
			
			if(qID.get(contig).containsKey(pos)){
				int id=qID.get(contig).get(pos);
				if(!removed.containsKey(id)){
					return id;
				}else
					return -1;
			}else{
				return -1;
			}
		}else{
			return -1;
		}	
	}
	
	public String get(int query){
		if(query<qIDList.size()&&query>=0&&!removed.containsKey(query)){
			Query q=qIDList.get(query);
			return contigsList.get(q.contig)+"_"+q.pos;	
		}else{
			return null;
		}
	}
	
	public ArrayList<Integer> add(ArrayList<String> queries){
		ArrayList<Integer> keys=new ArrayList<Integer>();
		for(int i=0;i<queries.size();i++){
			keys.add(add(queries.get(i)));
		}
		return keys;
	}
	
	public ArrayList<Integer> get(ArrayList<String> queries){
		ArrayList<Integer> keys=new ArrayList<Integer>();
		for(int i=0;i<queries.size();i++){
			int key=get(queries.get(i));
			if(key>-1)keys.add(key);
		}
		return keys;
	}
	
	public int add(String query){
		Query q=parse(query);
		int pos=q.pos;
		int contig=q.contig;
		if(qID.containsKey(contig)){
			if(qID.get(contig).containsKey(pos)){
				int id=qID.get(contig).get(pos);
				if(removed.containsKey(id)){
					removed.remove(id);
				}
				return id;
			}else{
				int id=qIDList.size();
				qID.get(contig).put(pos,id);
				qIDList.add(q);
				return id;
			}
		}else{
			HashMap<Integer,Integer> hm=new HashMap<Integer, Integer>();
			int id=qIDList.size();
			hm.put(pos, id);
			qID.put(contig,hm);
			qIDList.add(q);
			return id;
		}
	}
	
	private int addContigs(String contig){
		if(!contigs.containsKey(contig)){
			int id=contigsList.size();
			contigs.put(contig, id);
			contigsList.add(contig);
			return id;
		}else{
			return contigs.get(contig);
		}
	}
	
	public Query parse(String q){
		String[] split=q.split("_");
		String posS=split[split.length-1];
		int pos=Integer.parseInt(posS);
		String contigS=q.substring(0,q.length()-posS.length()-1);
		int contig=addContigs(contigS);
		return new Query(contig,pos);
		

	}
	
	class Query implements Serializable{
		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;
		int contig;
		int pos;
		public Query(int contig,int pos){
			this.contig=contig;
			this.pos=pos;
		}

	}
	 
}
