package GreedyExperimentalDesign;

import java.util.ArrayList;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;

public class GreedySearch {

	public GreedySearch(double[][] Xstd, DenseMatrix sinvmat, int[] indicT, int[] ending_indicT, Double[] objective_vals, String objective, int d0) {
//		System.out.println("GreedySearch: ready to begin " + d0);
		ObjectiveFunction objective_fun = null;
		if (objective.equals(GreedyExperimentalDesign.MAHAL)){
			objective_fun = new PropMahalObjective(sinvmat);
		} 
		else if (objective.equals(GreedyExperimentalDesign.ABS)){
			objective_fun = new AbsSumObjective();	
		}
		System.out.println("indicT: [" + Tools.StringJoin(indicT) + "]");
		beginSearch(Xstd, sinvmat, indicT, ending_indicT, objective_vals, d0, objective_fun);		
	}

	private void beginSearch(double[][] Xstd, DenseMatrix sinvmat, int[] indicT, int[] ending_indicT, Double[] objective_vals, int d0, ObjectiveFunction obj_fun) {
		int n = Xstd.length;
		int p = Xstd[0].length;		
		int nT = count(indicT, 1);
		
		Double obj_val = null;
		
		double min_obj_val = Double.MAX_VALUE;
		int iter = 0;
		while (true){
			iter++;
			System.out.println("GreedySearch: iter " + iter);			
			
			int[] indicTmin = null;
			int[] i_Ts = findIndicies(indicT, nT, 1);
			int[] i_Cs = findIndicies(indicT, n - nT, 0);
			for (int i_T : i_Ts){
				for (int i_C : i_Cs){
					int[] indicT_proposal = indicT.clone();
					//make the single switch
					indicT_proposal[i_T] = 0;
					indicT_proposal[i_C] = 1;
					
					ArrayList<double[]> XT = subsetMatrix(Xstd, nT, p, i_Ts, i_T, i_C); 
					ArrayList<double[]> XC = subsetMatrix(Xstd, nT, p, i_Cs, i_C, i_T); 

					obj_fun.setXTbar(colAvg(XT, p));
					obj_fun.setXCbar(colAvg(XC, p));
								
					obj_val = obj_fun.calc();
					
					if (obj_val < min_obj_val){
						indicTmin = indicT_proposal;
						indicT = indicT_proposal;
						min_obj_val = obj_val;
					}
				}				
			}
			//after searching through every possible switch, we didn't find anything, so break
			if (indicTmin == null){
				break;
			}			
		}		
		
		//search is over; ship back the data now
		for (int i = 0; i < indicT.length; i++){
			ending_indicT[i] = indicT[i];
		}
		objective_vals[d0] = obj_val;	
	}

	private DenseVector colAvg(ArrayList<double[]> X, int p) {
		int n = X.size();	
		double[] tally = new double[p];
		for (int i = 0; i < n; i++){
			for (int j = 0; j < p; j++){
				tally[j] += X.get(i)[j];
			}			
		}
		for (int j = 0; j < p; j++){
			tally[j] /= n;
		}
		return new DenseVector(tally);
	}

	private ArrayList<double[]> subsetMatrix(double[][] Xstd, int nT, int p, int[] indices, int i_remove, int i_add) {
		ArrayList<double[]> XstdT = new ArrayList<double[]>(nT);
		for (int i : indices){
			if (i != i_remove){
				XstdT.add(Xstd[i]);
			}			
		}
		XstdT.add(Xstd[i_add]);		
		return XstdT;
	}

	private int count(int[] indicT, int val) {
		int tally = 0;
		for (int i = 0; i < indicT.length; i++){
			if (indicT[i] == val){
				tally++;
			}			
		}
		return tally;
	}

	private int[] findIndicies(int[] vec, int n_val, int val) {
		int[] indicies = new int[n_val];
		int index = 0;
		for (int i = 0; i < vec.length; i++){
			if (vec[i] == val){
				System.out.println("index found at loc = " + i);
				indicies[index] = i;
				index++;
			}
		}
		return indicies;
	}
}
