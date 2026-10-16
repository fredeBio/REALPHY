package buildTreesFromSequenceData;

import java.io.Serializable;
import java.util.*;
import java.util.Map.Entry;


public class QueryBase implements Serializable{
		/**
	 * 
	 */

	private static final long serialVersionUID = 1L;
		//list of queryID columns filled in for each clash
		ArrayList<ArrayList<Integer>> queryIDs=new ArrayList<ArrayList<Integer>>();
		//same as for querIDs but the corresponding nucleotides instead of qIDs
		ArrayList<StringBuffer> baseColumns=new ArrayList<StringBuffer>();
		public void addAll(QueryBase qb){
			queryIDs.addAll(qb.queryIDs);
			baseColumns.addAll(qb.baseColumns);
		}
		public QueryBase(){
			
		}
		
		
		public QueryBase(ArrayList<ArrayList<Integer>> queryIDCols,ArrayList<StringBuffer> baseCols ){
			this.queryIDs=queryIDCols;
			this.baseColumns=baseCols;
		}
		/**
		 * Sorts the alignment by the position of the first reference. Recommended once all clashes are resolved.
		 */
		public void sort(){
			TreeMap<Integer,Integer> tm=new TreeMap<Integer, Integer>();
			for(int i=0;i<queryIDs.size();i++){
				tm.put(queryIDs.get(i).get(0),i);
			}
			ArrayList<ArrayList<Integer>> sortedqIDs=new ArrayList<ArrayList<Integer>>();
			ArrayList<StringBuffer> sortedBaseColumns=new ArrayList<StringBuffer>();
			Iterator<Entry<Integer,Integer>> it=tm.entrySet().iterator();
			while(it.hasNext()){
				Entry<Integer,Integer> e=it.next();
				int i=e.getValue();
				sortedqIDs.add(queryIDs.get(i));
				sortedBaseColumns.add(baseColumns.get(i));
				
			}
			queryIDs=sortedqIDs;
			baseColumns=sortedBaseColumns;
		}
		
		/**TODO
		 * Once the alignment is sorted gaps are added to the positions that are not present.
		 */
		public void addGaps(){
			int last=-1;
			int number=queryIDs.get(0).size();
			StringBuffer gaps=getGaps(number);
			ArrayList<ArrayList<Integer>> gappedqIDs=new ArrayList<ArrayList<Integer>>();
			ArrayList<StringBuffer> gappedBases=new ArrayList<StringBuffer>();
			for(int i=0;i<queryIDs.size();i++){
				int current=queryIDs.get(i).get(0);
				for(int j=last;j<current-1;j++){
					gappedqIDs.add(getQIDs(j,number));
					gappedBases.add(gaps);
				}
				gappedqIDs.add(queryIDs.get(i));
				gappedBases.add(baseColumns.get(i));
				last=current;
			}
			queryIDs=gappedqIDs;
			baseColumns=gappedBases;
		}
		
		private ArrayList<Integer> getQIDs(int pos,int number){
			ArrayList<Integer> list=new ArrayList<Integer>();
			for(int i=0;i<number;i++){
				list.add(pos);
			}
			return list;
		}
		
		private StringBuffer getGaps(int number){
			StringBuffer sb=new StringBuffer();
			for(int i=0;i<number;i++){
				sb.append('-');
			}
			return sb;
		}
		
		public QueryBase(ArrayList<Integer> queryIDColumn,StringBuffer bases ){
			this.queryIDs.add(queryIDColumn);
			this.baseColumns.add(bases);
		}
		public int length(){
			return queryIDs.size();
		}
		
	
}
