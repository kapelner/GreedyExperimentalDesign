package GreedyExperimentalDesign;

import no.uib.cipr.matrix.DenseMatrix;

public class GreedySearch {

	public GreedySearch(double[][] Xstd, DenseMatrix sinvmat, int[] indicT, int[] ending_indicT, Double[] objective_vals, String objective, int d0) {
		if (objective.equals(GreedyExperimentalDesign.MAHAL)){
			searchMahal(Xstd, sinvmat, indicT, ending_indicT, objective_vals, d0);		
		} 
		else if (objective.equals(GreedyExperimentalDesign.ABS)){
			searchAbs(Xstd, indicT, ending_indicT, objective_vals, d0);		
		}
	}

	private void searchAbs(double[][] xstd, int[] indicT, int[] ending_indicT, Double[] objective_vals, int d0) {

		
		//ship back the data now
//		for (int i = 0; i < indicT.length; i++){
//			ending_indicT[i] = indicT[i];
//		}
//		objective_vals[d0]		
	}

	private void searchMahal(double[][] xstd, DenseMatrix sinvmat, int[] indicT, int[] ending_indicT, Double[] objective_vals, int d0) {
		
		
		
		//ship back the data now
//		for (int i = 0; i < indicT.length; i++){
//			ending_indicT[i] = indicT[i];
//		}
//		objective_vals[d0]
	}



}
