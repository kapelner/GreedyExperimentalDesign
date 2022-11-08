package PairwiseMatchingExperimentalDesign;


import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;

import ExperimentalDesign.MultipleSearchExperimentalDesigns;

public class PairwiseMatchingExperimentalDesign extends MultipleSearchExperimentalDesigns {

	private int[][] binary_pairs;

	private HashMap<Integer, int[]> unique_allocations;
	
	public void beginSearch(){
		super.beginSearch();
		
		unique_allocations = new HashMap<Integer, int[]>(max_designs);
		
		for (int d = 0; d < max_designs; d++){
	    	search_thread_pool.execute(new Runnable(){
				public void run() {
					if (search_stopped.get()) {
						return;
					}
					while (true) {
						int[] w = generateAllocation();
						Integer h = Arrays.hashCode(w);
						if (unique_allocations.get(h) == null) {
							unique_allocations.put(h, w);
							break;
						}
					}
					
					num_completed.getAndIncrement();
//					System.out.println("did one num_completed: " + num_completed.get());
					
				}
			});
		}	
		afterBeginSearch();
	}
	
	protected void afterBeginSearch() {
		super.afterBeginSearch();
		
//		UniqueAllocation[] unique_allocations_array = (UniqueAllocation[])unique_allocations.toArray();
//		unique_allocations = null; //why not cleanup right away?
//		
		//now copy into the allocations into the canonical form
		Object[] allocations = unique_allocations.values().toArray();
		for (int d = 0; d < max_designs; d++) {
			ending_indicTs[d] = (int[])allocations[d];
		}
	}
	
	protected int[] generateAllocation() {
		int[] w = new int[n];
		for (int m = 0; m < (n / 2); m++) {
			int[] indicies_pair = binary_pairs[m];
			if (rand_obj.nextBoolean()) {
				w[indicies_pair[0]] = 1;
				w[indicies_pair[1]] = 0;
			} else {
				w[indicies_pair[1]] = 1;
				w[indicies_pair[0]] = 0;				
			}
		}
		return w;
	}

	public void setMatchPairIndicies(int i, int[] pair) {
		if (binary_pairs == null) {
			binary_pairs = new int[n / 2][];
		}		
		binary_pairs[i] = pair;
	}	
}
