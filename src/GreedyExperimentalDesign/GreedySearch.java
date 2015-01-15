package GreedyExperimentalDesign;

import java.util.ArrayList;

import no.uib.cipr.matrix.DenseMatrix;

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
		int nT = Tools.count(indicT, 1);
		System.out.println("beginSearch: nT = " + nT + " and nC = " + (n - nT));
		
		Double obj_val = null;
		
		double min_obj_val = Double.MAX_VALUE;
		int iter = 0;
		while (true){
			iter++;
						
			
			int[] indicTmin = null;
			int[] i_Ts = Tools.findIndicies(indicT, nT, 1);
			int[] i_Cs = Tools.findIndicies(indicT, n - nT, 0);
			
			System.out.println("GreedySearch: iter " + iter + " #i_Ts: " + i_Ts.length + " #i_Cs: " + i_Cs.length);
			for (int i_T : i_Ts){
				for (int i_C : i_Cs){
					
					int[] indicT_proposal = indicT.clone();
					//make the single switch
					indicT_proposal[i_T] = 0;
					indicT_proposal[i_C] = 1;
					
					ArrayList<double[]> XT = Tools.subsetMatrix(Xstd, nT, p, i_Ts, i_T, i_C); 
					ArrayList<double[]> XC = Tools.subsetMatrix(Xstd, nT, p, i_Cs, i_C, i_T); 

					obj_fun.setXTbar(Tools.colAvg(XT, p));
					obj_fun.setXCbar(Tools.colAvg(XC, p));
								
					obj_val = obj_fun.calc();
					
					if (obj_val < min_obj_val){
						indicTmin = indicT_proposal;
						indicT = indicT_proposal;
						min_obj_val = obj_val;
						System.out.println("switched i_T " + i_T + " and i_C " + i_C);
						System.out.println("min_obj_val " + min_obj_val);
					}
				}				
			}
			System.out.println("end of double loop");
			System.out.println("indicT: " + indicT + " indicTmin: " + indicTmin);
			//after searching through every possible switch, we didn't find anything, so break
			if (indicTmin == null){
//				System.out.println("break");
				break;
			}
			System.out.println("after while true");
		}		
		
		//search is over; ship back the data now
		for (int i = 0; i < indicT.length; i++){
			ending_indicT[i] = indicT[i];
		}
		objective_vals[d0] = obj_val;
		System.out.println("SEARCH DONE obj_val " + objective_vals[d0]);
	}

}
