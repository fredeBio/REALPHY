package buildTreesFromSequenceData;

import java.io.File;
import java.util.HashMap;

import solexa.PointSubstitutions;




public interface GetPolymorphisms {

      HashMap<String /*geneID*/,HashMap<Integer/*position*/,Polymorph/*poly*/>> adjustGP(PointSubstitutions pss,boolean subInfo,String name);

	  void checkCoverageGP(PointSubstitutions pss,String strain);
	public Clashes calculateColumns();


	
	public File writeSequences();
	

}
