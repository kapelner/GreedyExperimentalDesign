package KarpExperimentalDesign;

import java.util.Collections;
import java.util.Comparator;


public class KarpDesignSearcherBalanced extends KarpDesignSearcher {
	
	public KarpDesignSearcherBalanced(double[][] Xstd) {
		super(Xstd);
	}
	
	
	protected class ObsBundleCompareBalanced implements Comparator<ObsBundle> {
		@Override
		public int compare(ObsBundle o1, ObsBundle o2) {
			if (o1.size() < o2.size()){
				return -1;
			}
			else if (o1.size() > o2.size()){
				return 1;
			}
			else {
				if (o1.x_val < o2.x_val){
					return 1;
				}
				else if (o1.x_val > o2.x_val){
					return -1;
				}
				else {
					return 0;
				}				
			}
		}		
	}
	
	public void sortObsBundles(){
		Collections.sort(obs_bundles, new ObsBundleCompareBalanced());
	}
}
