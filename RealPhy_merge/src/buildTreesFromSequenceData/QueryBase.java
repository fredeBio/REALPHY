package buildTreesFromSequenceData;

import java.io.Serializable;
import java.util.*;


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
		
		public QueryBase(ArrayList<Integer> queryIDColumn,StringBuffer bases ){
			this.queryIDs.add(queryIDColumn);
			this.baseColumns.add(bases);
		}
		public int length(){
			return queryIDs.size();
		}
	
}
