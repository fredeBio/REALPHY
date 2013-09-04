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

	/*public class Pair{
		int count;
		ArrayList<Integer> refs=new ArrayList<Integer>();
		public Pair(int ref){
			count=1;
			refs.add(ref);
		}
		private void increaseCount(){
			this.count++;

		}
		public Pair addRef(int ref){
			this.refs.add(ref);
			increaseCount();
			return this;
		}
		
		
	}*/
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	/**
	 * is filled out after program has found all the clashes, contains the position of each queryID in the column list, when implementing reference position only a single positional reference is possible
	 * this means that a column ID can only refer to a single position in the query ID array and this is the position where this position was the reference
	 */

	HashMap</*position in column*/Integer,HashMap<Integer/*position in genome*/,ArrayList<Integer>>> posCol=new HashMap<Integer,HashMap<Integer, ArrayList<Integer>>>();
	//filled out after program has found all clashes, updated every time a clash has been resolved, indicates the positions in the AL that have been resolved
	HashMap<Integer,Boolean> resolved=new HashMap<Integer, Boolean>();
	//number of positions of unresolved clashes, is filled in after all clashes have been found and gradually decreases as clashes get resolved
	ArrayList<Integer> unresolved=new ArrayList<Integer>();
	
	//keeps track of the queryIDs that are currently found in the clash class
	//HashMap<String,Boolean> clashes=new HashMap<String, Boolean>();
	//class to combine list of queryIds and list of base columns
	QueryBase qb=new QueryBase();
	
	//counts the number of identical columns
	HashMap<ArrayList<Integer>,ArrayList<Integer>> count=new HashMap<ArrayList<Integer>, ArrayList<Integer>>();
	
	//needed for alignment output!
	ArrayList<String> idents;
	
	public Clashes setIdents(ArrayList<String> idents){
		this.idents=idents;
		return this;
	}
	public Clashes(){
		
	}
	
	public ArrayList<Integer> getQIDColumn(int pos){
		return qb.queryIDs.get(pos);
	}
	
	public StringBuffer getBaseColumn(int pos){
		return qb.baseColumns.get(pos);
	}
	
	public void addColumn(ArrayList<Integer> queryIDColumn,StringBuffer bases,int ref){

		if(!count.containsKey(queryIDColumn)){
			ArrayList<Integer> p=new ArrayList<Integer>();
			p.add(ref);
			count.put((queryIDColumn), p);
			qb.queryIDs.add(queryIDColumn);
			qb.baseColumns.add(bases);
			unresolved.add(qb.queryIDs.size()-1);;

		}else{

			count.get(queryIDColumn).add(ref);
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
			if(unresolved.size()%10000==0){
				System.out.println("Sites still to resolve: "+ unresolved.size());
			}
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
		ArrayList<Integer> qIDs=qb.queryIDs.get(rand);
		ArrayList<StringBuffer> basesCons=new ArrayList<StringBuffer>();
		ArrayList<ArrayList<Integer>> qIDCons=new ArrayList<ArrayList<Integer>>();
		for(int i=0;i<qIDs.size();i++){		
			if(posCol.get(i).containsKey(qIDs.get(i))){
				ArrayList<Integer> pos=posCol.get(i).get(qIDs.get(i));
				for(int j=0;j<pos.size();j++){
					if(!resolved.containsKey(pos.get(j))){

						ArrayList<Integer> columns=qb.queryIDs.get(pos.get(j));
						StringBuffer bases=qb.baseColumns.get(pos.get(j));
						ArrayList<Integer> p=count.get((columns));
						for(int k=0;k<p.size();k++){
							if(p.get(k)==i){
								basesCons.add(bases);
								qIDCons.add(columns);
								p.remove(k);
								break;
							}

						}
						if(p.size()==0){
							resolved.put(pos.get(j),true);;
						}

					}

				}
			}
		}
		QueryBase consense=null;
		if(basesCons.size()>0){
//			int size=basesCons.size();
			
			consense=getConsensus(basesCons,qIDCons);
//			if(consHM.containsKey(consense.queryIDs.get(0).get(0))){
//				System.err.println("PROBLEM "+consense.queryIDs.get(0).get(0));
//				//System.exit(-1);
//
//				consHM.put(consense.queryIDs.get(0).get(0), consHM.get(consense.queryIDs.get(0).get(0))+1);
//			}else{
//				consHM.put(consense.queryIDs.get(0).get(0),1);
//			}
			/*if(basesCons.size()>1&&!basesCons.get(0).equals(basesCons.get(1))){
				System.out.println("all: "+qIDCons);
				System.out.println("all: "+basesCons);
				System.out.println("consense: "+consense.baseColumns);
			}*/
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
			ArrayList<Integer> col=newClashes.getQIDColumn(i);
			ArrayList<Integer> p=newClashes.count.get((col));
			StringBuffer baseColumn=newClashes.getBaseColumn(i);
			for(int j=0;j<p.size();j++){
				addColumn(col,baseColumn ,p.get(j));
				
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
	
	private boolean isEqual(ArrayList<StringBuffer> bases){
		String ref=bases.get(0).toString().toUpperCase();
		for(int i=1;i<bases.size();i++){
			String test=bases.get(i).toString().toUpperCase();
			if(!ref.equals(test)){
				return false;
			}
		}
		return true;
	}
	
	private QueryBase getConsensus(ArrayList<StringBuffer> bases,ArrayList<ArrayList<Integer>> qIDs){
		ArrayList<Integer> qIDsCons=new ArrayList<Integer>();
		StringBuffer basesCons=new StringBuffer();
		HashMap<Integer,Integer> colHash=new HashMap<Integer, Integer>();
		if(!isEqual(bases)){
			
			ArrayList<String> baseCols=initbaseCols(bases);
			int length=baseCols.get(0).length();
			for(int j=0;j<length;j++){
				HashMap<Character,Integer> baseHash=new HashMap<Character, Integer>();
				for(int i=0;i<qIDs.size();i++){
					if(j==0){
						System.err.println(qIDs.get(i));
						System.err.println(baseCols.get(i));
					}
					String baseCol=baseCols.get(i);
					ArrayList<Integer> col=qIDs.get(i);
					char c=baseCol.charAt(j);

					enterItem(baseHash,c);
					if(j<col.size()){
						Integer qID=col.get(j);
						enterItem(colHash,qID);
					}
				}
				qIDsCons.add(getMax(colHash));
				basesCons.append(getMax(baseHash));
			}
			System.err.println("next");

			return new QueryBase(qIDsCons,basesCons);
		}else{
			return new QueryBase(qIDs.get(0),bases.get(0));
		}
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
	 * sets the position of the reference queryIDs from a given starting position in the qb alignment, now it is possible to determine the reference column for a given query mapping
	 * @param start
	 */
	private void createPosCol(int start){
		for(int i=start;i<qb.queryIDs.size();i++){
			if(i%10000==0){
				System.out.println("Created positional index for "+i+" out of "+qb.queryIDs.size());
			}
			ArrayList<Integer> col=qb.queryIDs.get(i);
			ArrayList<Integer> refs=count.get(col);
			for(int j=0;j<refs.size();j++){
				int ref=refs.get(j);

				int colPos=ref;
				int genomePos=col.get(ref);

				if(posCol.containsKey(colPos)){
					if(posCol.get(colPos).containsKey(genomePos)){
						//System.err.println("There is something wrong. It is not possible that one query position refers to two different reference positions!");
						//System.exit(-1);
						//this is possible since the reference is also mapped to the reference...
						posCol.get(colPos).get(genomePos).add(i);
					}else{
						ArrayList<Integer> temp=new ArrayList<Integer>();
						temp.add(i);
						posCol.get(colPos).put(genomePos,temp);
					}
				}else{
					ArrayList<Integer> temp=new ArrayList<Integer>();
					temp.add(i);
					HashMap<Integer,ArrayList<Integer>> hm=new HashMap<Integer, ArrayList<Integer>>();
					hm.put(genomePos,temp);
					posCol.put(colPos,hm);
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
		resolved.sort();

		//resolved.addGaps();
		for(int i=0;i<resolved.length();i++){
			alg.addColumn(resolved.baseColumns.get(i).toString());
		}
		Fasta.write(alg.toFasta(),out);
		
		//print();
	}
	
	public Clashes(ArrayList<File> clashFiles){

		mergeFiles(clashFiles);
	}
	
	/**
	 * automatically adds new set of columns stored in a list of files to current data structure
	 * @param clashObjectFiles
	 */
	public  void mergeFiles(ArrayList<File> clashObjectFiles){
		int size=clashObjectFiles.size();
		for(int i=0;i<size;i++){
			System.out.println(Runtime.getRuntime().totalMemory());
			System.out.println("Loading and merging "+clashObjectFiles.get(i).getName()+" "+i+" out of "+size+".");
			System.out.println("Current clash size "+this.size());
			
			Clashes c=(Clashes)ObjectIO.readObject(clashObjectFiles.get(i));
			System.out.println("Adding reference: "+c.count.get(c.getQIDColumn(0)).get(0));
			//printColumn(c.getQIDColumn(0));
			if(this.idents!=null){
				if(isEqual(c.idents,this.idents)){


				}else{
					System.out.println("Idents are not equal!");
					printIdents(c.idents);

				}
			}else{
				this.idents=c.idents;
			}
			
			this.addAll(c);
			
		}
	}
	
	private void printColumn(ArrayList<Integer> id){
		for(int i=0;i<id.size();i++){
			System.out.println(id.get(i));
		}
	}
	
	private void printIdents(ArrayList<String> id){
		for(int i=0;i<id.size();i++){
			System.out.println(id.get(i));
		}
	}
	
	private boolean isEqual(ArrayList<String> id1,ArrayList<String> id2){
		if(id1.size()==id2.size()){
			for(int i=0;i<id1.size();i++){
				if(!id1.get(i).equals(id2.get(i))){
					return false;
				}
			}
		}else{
			return false;
		}
		return true;
	}
}
