package buildTreesFromSequenceData;

import java.io.File;
import java.io.Serializable;
import java.util.*;
import java.util.Map.Entry;

import pairwiseAlignment.NeedlemanWunsch;

import util.DNAmanipulations;
import util.Fasta;
import util.ObjectIO;
import util.phylogenetics.Alignment;

public class Clashes implements Serializable{


	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	//is filled out after program has found all the clashes, contains the position of each queryID in the column list
	HashMap<String,ArrayList<Integer>> posCol=new HashMap<String, ArrayList<Integer>>();
	//filled out after program has found all clashes, updated every time a clash has been resolved, indicates the positions in the AL that have been resolved
	HashMap<Integer,Boolean> resolved=new HashMap<Integer, Boolean>();
	//number of positions of unresolved clashes, is filled in after all clashes have been found and gradually decreases as clashes get resolved
	ArrayList<Integer> unresolved=new ArrayList<Integer>();
	
	//keeps track of the queryIDs that are currently found in the clash class
	//HashMap<String,Boolean> clashes=new HashMap<String, Boolean>();
	//class to combine list of queryIds and list of base columns
	QueryBase qb=new QueryBase();
	
	//counts the number of identical columns
	HashMap<ArrayList<String>,Integer> count=new HashMap<ArrayList<String>, Integer>();
	
	//needed for alignment output!
	ArrayList<String> idents=new ArrayList<String>();
	
	public Clashes(ArrayList<String> idents){
		this.idents=idents;
	}
	
	public ArrayList<String> getQIDColumn(int pos){
		return qb.queryIDs.get(pos);
	}
	
	public StringBuffer getBaseColumn(int pos){
		return qb.baseColumns.get(pos);
	}
	
	public void addColumn(ArrayList<String> qIDCol,StringBuffer bases){

		if(!count.containsKey(qIDCol)){
			count.put((qIDCol), 1);
			qb.queryIDs.add(qIDCol);
			qb.baseColumns.add(bases);
			unresolved.add(qb.queryIDs.size()-1);;

//			for(int i=0;i<qIDCol.size();i++){
//				clashes.put(qIDCol.get(i),true);
//			}
		}else{
			//String q=toString(qIDCol);
			count.put(qIDCol, count.get(qIDCol)+1);
		}
	}
	
//	public boolean isClash(ArrayList<String> qIDs){
//		for(int i=0;i<qIDs.size();i++){
//			if(clashes.containsKey(qIDs.get(i))) return true;
//		}
//		return false;
//	}
	
	public QueryBase resolveClashes(){
		createPosCol();
		QueryBase resolved=new QueryBase();
		while(unresolved.size()>0){
			int rand=(int)(Math.random()*unresolved.size());
			QueryBase result=resolveClash(unresolved.get(rand));
			if(result!=null){
				resolved.addAll(result);
			}
			unresolved.remove(rand);
		}
		return resolved;
	}
	
//	public void addColumns(QueryBase newQB){
//		for(int i=0;i<newQB.queryIDs.size();i++){
//			addColumn(newQB.queryIDs.get(i),newQB.baseColumns.get(i));
//		}
//	}
	
//	private void print(ArrayList<String> col){
//		for(int i=0;i<col.size();i++){
//			System.out.println(col.get(i));
//		}
//	}
	
	//takes a random column out of the list of clashed columns as well as all other columns in which one of the involved column IDs is found
	//it then calculates the consensus of all these columns and stores it in Column format
	private QueryBase resolveClash(int rand){
		ArrayList<String> qIDs=qb.queryIDs.get(rand);
		ArrayList<StringBuffer> basesCons=new ArrayList<StringBuffer>();
		ArrayList<ArrayList<String>> qIDCons=new ArrayList<ArrayList<String>>();
		for(int i=0;i<qIDs.size();i++){			
			ArrayList<Integer> positions=posCol.get(qIDs.get(i));
			for(int j=0;j<positions.size();j++){
				int pos=positions.get(j);
				if(!resolved.containsKey(pos)){
					resolved.put(pos,true);
					ArrayList<String> columns=qb.queryIDs.get(pos);
					StringBuffer bases=qb.baseColumns.get(pos);
					int c=count.get((columns));
					for(int k=0;k<c;k++){
						basesCons.add(bases);
						qIDCons.add(columns);
					}

				}
			}
		}
		QueryBase consense=null;
		if(basesCons.size()>0){
			
			consense=getConsensus(basesCons,qIDCons);
			if(basesCons.size()>1&&!basesCons.get(0).equals(basesCons.get(1))){
				System.out.println("all: "+qIDCons);
				System.out.println("all: "+basesCons);
				System.out.println("consense: "+consense.baseColumns);
			}
		}else{
			/*		System.out.println(rand);
			System.out.println(qb.queryIDs.get(rand));
			for(int i=0;i<qIDs.size();i++){			
				ArrayList<Integer> positions=posCol.get(qIDs.get(i));
				for(int j=0;j<positions.size();j++){
					int pos=positions.get(j);
					System.out.println(qb.queryIDs.get(pos));
					System.out.println(qb.baseColumns.get(pos));
				}
			}
							System.exit(-1);
*/
		}
		return consense;
	}
	
//	private String toString(ArrayList<String> l){
//		StringBuffer sb=new StringBuffer();
//		for(int i=0;i<l.size();i++){
//			sb.append(l.get(i));
//		}
//		return sb.toString();
//	}
	
	/**
	 * adds more clashes, resets all "resolved" clashes
	 * @param newClashes
	 */
	public void addAll(Clashes newClashes){
		resolved=new HashMap<Integer, Boolean>();

		for(int i=0;i<newClashes.size();i++){
			ArrayList<String> col=newClashes.getQIDColumn(i);
			int c=newClashes.count.get((col));
			for(int j=0;j<c;j++){
				addColumn(newClashes.getQIDColumn(i), newClashes.getBaseColumn(i));
				
			}
		}
		//System.out.println("Clash Size: "+this.size()+" "+unresolved.size());

	}
	
	/**
	 * prints the number of columns that are currently clashing
	 * @return
	 */
	public int size(){
		return qb.length();
	}
	
	private ArrayList<String> initbaseCols(ArrayList<StringBuffer> al){
		ArrayList<String> newAl=new ArrayList<String>();
		String ref;
		if(al.size()>0){
			ref=al.get(0).toString().toUpperCase();
			newAl.add(ref);
			for(int i=1;i<al.size();i++){
				String s=al.get(i).toString().toUpperCase();
				if(!s.equals(ref)){
					String srev=DNAmanipulations.complement(s);
					double idN= NeedlemanWunsch.getPairwiseIdentity(s, ref);
					double idRev=NeedlemanWunsch.getPairwiseIdentity(srev, ref);
					if(idN>idRev){
						newAl.add(s);
					}else{
						newAl.add(srev);
					}
				}else{
					newAl.add(s);
				}
			}
		}
		return newAl;
	}
	
	/**
	 * calculates a consensus for columns and queryIDs, assumes that it receives a sorted list!
	 * @param bases
	 * @param qIDs
	 * @return
	 */
	
	private QueryBase getConsensus(ArrayList<StringBuffer> bases,ArrayList<ArrayList<String>> qIDs){
		ArrayList<String> qIDsCons=new ArrayList<String>();
		StringBuffer basesCons=new StringBuffer();
		HashMap<String,Integer> colHash=new HashMap<String, Integer>();
		ArrayList<String> baseCols=initbaseCols(bases);

			for(int j=0;j<baseCols.get(0).length();j++){
				HashMap<Character,Integer> baseHash=new HashMap<Character, Integer>();
				for(int i=0;i<qIDs.size();i++){
					String baseCol=baseCols.get(i);
					ArrayList<String> col=qIDs.get(i);
					char c=baseCol.charAt(j);
					
					enterItem(baseHash,c);
					if(j<col.size()){
						String qID=col.get(j);
						enterItem(colHash,qID);
					}
				}
				qIDsCons.add(getMax(colHash));
				basesCons.append(getMax(baseHash));
		}
		return new QueryBase(qIDsCons,basesCons);
	}
	
	private <T> void enterItem(HashMap<T,Integer> hm,T item){
		if(hm.containsKey(item)){
			hm.put(item, hm.get(item)+1);
		}else{
			hm.put(item, 1);
		}
	}
	
	private <T> T getMax(HashMap< T,Integer> hm){
		Iterator<Entry<T,Integer>> it=hm.entrySet().iterator();
		int max=Integer.MIN_VALUE;
		T maxO=null;
		while(it.hasNext()){
			Entry<T,Integer> e=it.next();
			if(max<e.getValue()){
				max=e.getValue();
				maxO=e.getKey();
			}
		}
		return maxO;
	}
	
	/**
	 * sets the position of all queryIDs from a given starting position in the qb alignment
	 * @param start
	 */
	private void createPosCol(int start){
		for(int i=start;i<qb.queryIDs.size();i++){
			ArrayList<String> col=qb.queryIDs.get(i);
			for(int j=0;j<col.size();j++){
				if(posCol.containsKey(col.get(j))){
					posCol.get(col.get(j)).add(i);
				}else{
					ArrayList<Integer> al=new ArrayList<Integer>();
					al.add(i);
					posCol.put(col.get(j),al);
				}
			}
		}
		
	}
	private void createPosCol(){
		createPosCol(0);
		
	}
	
	/**
	 * Turns the clash into an alignment, AFTER resolving the clashes!!!
	 * **/

	public void printAlignment(File out){
		QueryBase resolved=resolveClashes();
		Alignment alg=new Alignment();
		alg.changeIdents(idents);
		for(int i=0;i<resolved.length();i++){
			alg.addColumn(resolved.baseColumns.get(i).toString());
		}
		Fasta.write(alg.toFasta(),out);
	}
	
	/**
	 * automatically adds new set of columns stored in a list of files to current data structure
	 * @param clashObjectFiles
	 */
	public  void mergeFiles(ArrayList<File> clashObjectFiles){
		for(int i=0;i<clashObjectFiles.size();i++){
			this.addAll((Clashes)ObjectIO.readObject(clashObjectFiles.get(i)));
			
		}
	}

}
