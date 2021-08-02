package GreedyExperimentalDesign;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;
import java.util.concurrent.atomic.AtomicBoolean;

import ExperimentalDesign.*;
import ObjectiveFunctions.*;

public class GreedySearch {

	private double[][] Xscaled;
	private double[][] X;
	private int nT;

	public GreedySearch(
		double[][] X, 
		double[][] Sinvmat, 
		HashMap<Integer, int[]> legal_pairs,
		double[][] Kgram,
		int[] indicT, 
		int[] ending_indicT, 
		ArrayList<int[]> switched_pairs,
		ArrayList<Double> min_obj_val_by_iteration,
		ArrayList<double[]> xbardiffjs_by_iteration,
		Double[] objective_vals, 
		Integer[] num_iters, 
		String objective, 
		int d0, 
		boolean semigreedy, //purely experimental... we didn't see any gain in this
		boolean diagnostics, 
		Integer max_iters, 
		Random r, 
		AtomicBoolean search_stopped, 
		HashMap<Integer, double[][]> Kgrams, 
		HashMap<Integer, Double> max_reduction_log_obj_vals, 
		double[] kernel_weights, 
		Double maximum_gain_scaling, 
		ArrayList<double[]> kernel_obj_values
	) {

		
//		System.out.println("GreedySearch: ready to begin " + d0);
		nT = Tools.count(indicT, 1);
		
		ObjectiveFunction obj_fun = null;
		int p = 0;
		int n = 0;
		if (objective.equals(ObjectiveFunction.KER)){
			obj_fun = new KernelObjective(Kgram);	
			n = Kgram.length;
			((KernelObjective)obj_fun).setW(indicT);	
		} else if (objective.equals(ObjectiveFunction.MUL_KER_PCT)) {
			obj_fun = new MultipleKernelObjectiveFunction(Kgrams, max_reduction_log_obj_vals, kernel_weights, maximum_gain_scaling, kernel_obj_values);	
			n = Kgrams.get(0).length; 
			((MultipleKernelObjectiveFunction)obj_fun).setW(indicT);
			((MultipleKernelObjectiveFunction)obj_fun).setInitialObjVals();
			if (diagnostics) {
				((MultipleKernelObjectiveFunction)obj_fun).resetKernelSum(); //waste, yes... but cleanest way to do it
				((MultipleKernelObjectiveFunction)obj_fun).calcKernelObjDiagnostics();
			}
			((MultipleKernelObjectiveFunction)obj_fun).resetKernelSum(); //waste, yes... but cleanest way to do it
			System.out.println("MultipleKernelObjectiveFunction GreedySearch #" + d0 + " ready to begin");
		} else {
			this.X = X;
			n = X.length;
			p = X[0].length;		
			
			createScaledXstd();
			
			if (objective.equals(ObjectiveFunction.MAHAL)){
				obj_fun = new MahalObjective(Sinvmat, n);
			} 
			else if (diagnostics && objective.equals(ObjectiveFunction.ABS)){
				obj_fun = new AbsSumObjectiveWithDiagnostics();	
			}
			else if (objective.equals(ObjectiveFunction.ABS)){
				obj_fun = new AbsSumObjective();
			}
		}
//		System.out.println("beginSearch d0:" + (d0 + 1) + " nT = " + nT + " and nC = " + (n - nT));
		
//		int[] i_Tss = Tools.findIndicies(indicT, nT, 1);
//		System.out.println("i_Ts " + Tools.StringJoin(i_Tss));
		
		Double obj_val = null;		
		
		
		int[] i_Ts = Tools.findIndicies(indicT, nT, 1);
//		System.out.println("i_Ts " + Tools.StringJoin(i_Ts));
		int[] i_Cs = Tools.findIndicies(indicT, n - nT, 0);
		
		ArrayList<double[]> XT = null;
		ArrayList<double[]> XC = null;
		double[] avg_Ts = null;
		double[] avg_Cs = null;
		
		if (obj_fun instanceof SimpleAverageObjectiveFunction) {
			XT = Tools.subsetMatrix(X, i_Ts); 
			XC = Tools.subsetMatrix(X, i_Cs);
			avg_Ts = Tools.colAvg(XT, p);
			avg_Cs = Tools.colAvg(XC, p);			
			((SimpleAverageObjectiveFunction)obj_fun).setXTbar(avg_Ts);
			((SimpleAverageObjectiveFunction)obj_fun).setXCbar(avg_Cs);
		}
		double min_obj_val = obj_fun.calc(false); //start at wherever we begin
			
		if (diagnostics){
//			System.out.println("calculating objective function for first time");
			min_obj_val_by_iteration.add(min_obj_val);
//			System.out.println("  iter 0 obj_val = " + min_obj_val);
		}
		
		int iter = 0;
		while (true){
			if (objective.equals(ObjectiveFunction.MUL_KER_PCT)) {
				System.out.println("    iter " + iter);
			}
			iter++;
//			System.out.println("iter++ " + iter);
			
			int[] indicTmin = null;
			double[] xbardiffjs = null;
//			System.out.println("indicTmin " + indicTmin);
			
			
			
			i_Ts = Tools.findIndicies(indicT, nT, 1);
//			System.out.println("i_Ts " + Tools.StringJoin(i_Ts));
			i_Cs = Tools.findIndicies(indicT, n - nT, 0);
//			System.out.println("i_Cs " + Tools.StringJoin(i_Cs));
			if (semigreedy){ //gotta randomize otherwise inefficient
				i_Ts = Tools.fisherYatesShuffle(i_Ts, r);
				i_Cs = Tools.fisherYatesShuffle(i_Cs, r);
			}
			//build the first avg vectors for speed

			if (obj_fun instanceof SimpleAverageObjectiveFunction) {
				XT = Tools.subsetMatrix(X, i_Ts); 
				XC = Tools.subsetMatrix(X, i_Cs); 
	
				avg_Ts = Tools.colAvg(XT, p);
				avg_Cs = Tools.colAvg(XC, p);
			}
//			System.out.println("INIT XTbar: " + Tools.StringJoin(avg_Ts, ","));
//			System.out.println("INIT XCbar: " + Tools.StringJoin(avg_Cs, ","));
			int[] switched_pair = new int[2];
			

						
//			System.out.println("iter " + iter + " #i_Ts: " + i_Ts.length + " #i_Cs: " + i_Cs.length);
			HashMap<Integer, int[]> Ts_to_Cs = setupIndicies(legal_pairs, i_Ts, i_Cs);
			
			indices_loop: {
				for (int i_T : Ts_to_Cs.keySet()){
					for (int i_C : Ts_to_Cs.get(i_T)){
//						System.out.println("   i_T " + i_T + " i_C " + i_C);
						
						

						
						int[] indicT_proposal = indicT.clone();
						//make the single switch
						indicT_proposal[i_T] = 0; //i_T is the new control
						indicT_proposal[i_C] = 1; //i_C is the new treatment
						
						if (objective.equals(ObjectiveFunction.KER)){	
							((KernelObjective)obj_fun).setSwitch(i_T, i_C);						
						} else if (objective.equals(ObjectiveFunction.MUL_KER_PCT)) {
							((MultipleKernelObjectiveFunction)obj_fun).setSwitch(i_T, i_C);		
						} else {
//							System.out.println("   updateAvgVec");
							updateAvgVec(avg_Ts, i_T, i_C, nT);
							((SimpleAverageObjectiveFunction)obj_fun).setXTbar(avg_Ts);
							
							updateAvgVec(avg_Cs, i_C, i_T, n - nT);
							((SimpleAverageObjectiveFunction)obj_fun).setXCbar(avg_Cs);
//							System.out.println("set XTbar and XCbar");
						}

						
						//calculate our objective function (according to the user's specification)
//						System.out.println("calculating objective function for iter " + iter + " i_T = " + i_T + " i_C = " + i_C);
						obj_val = obj_fun.calc(false);
						

//						System.out.println("  i_T = " + i_T + " i_C = " + i_C + " obj_val = " + obj_val);
						
						if (obj_val < min_obj_val){
							indicTmin = indicT_proposal;
//							System.out.println("best indicT so far " + Tools.StringJoin(indicTmin));
							min_obj_val = obj_val;
//							System.out.println("switched i_T " + i_T + " and i_C " + i_C);							
//							System.out.println("min_obj_val " + min_obj_val + " for iter " + iter);
							
							


							if (diagnostics){
								switched_pair[0] = i_T;
								switched_pair[1] = i_C;
								if (objective.equals(ObjectiveFunction.ABS)){
									xbardiffjs = ((AbsSumObjectiveWithDiagnostics)obj_fun).getXbardiffjs().clone();
								}
							}
							
							if (semigreedy){ //semigreedy means as soon as we find improvement, we ditch
								break indices_loop;
							}							
						}
						
						//reset the avg vecs
						if (obj_fun instanceof SimpleAverageObjectiveFunction){
							updateAvgVec(avg_Ts, i_C, i_T, nT);
							updateAvgVec(avg_Cs, i_T, i_C, n - nT);	
							
						}
					}	
				}
			}
//			System.out.println("after indices loop");
			
			//we've finished one iteration by checking every possible switch
			//record this switch only if it is a real switch
			if (diagnostics && indicTmin != null){
				switched_pairs.add(switched_pair);
				min_obj_val_by_iteration.add(min_obj_val);
				if (objective.equals(ObjectiveFunction.ABS)){
					xbardiffjs_by_iteration.add(xbardiffjs);
				}
			}
			
			//after searching through every possible switch, we didn't find anything, so break
			if (indicTmin == null){
				break;
			}
			//otherwise - continue and update the binary search vector
			else {
				indicT = indicTmin;
			}
			
			if (objective.equals(ObjectiveFunction.KER)){	
				((KernelObjective)obj_fun).resetKernelSum();
				((KernelObjective)obj_fun).setW(indicT);
				obj_fun.calc(false); //caches the current objective value
//				((KernelObjective)obj_fun).setPermanentSwitch(switched_pair[0], switched_pair[1]);						
			} else if (objective.equals(ObjectiveFunction.MUL_KER_PCT)){	
				if (diagnostics) {
					((MultipleKernelObjectiveFunction)obj_fun).calcKernelObjDiagnostics();
				}
				((MultipleKernelObjectiveFunction)obj_fun).resetKernelSum();
				((MultipleKernelObjectiveFunction)obj_fun).setW(indicT);
				obj_fun.calc(false); //caches the current objective value
//				((KernelObjective)obj_fun).setPermanentSwitch(switched_pair[0], switched_pair[1]);						
			}
//			
//			System.out.println("  iter " + iter + " obj_val = " + min_obj_val);

			//we can also be done if we hit our upper limit of iterations
			if (max_iters != null && max_iters == iter){
				break;
			}
//			System.out.println("before search_stopped");
			if (search_stopped.get()) {
				break;
			}
		}	
//		System.out.println("after while true");
		
		//search is over; ship back the data now
		for (int i = 0; i < indicT.length; i++){
			ending_indicT[i] = indicT[i];
		}
//		System.out.println("ending_indicT " + Tools.StringJoin(ending_indicT));
		objective_vals[d0] = min_obj_val;
		num_iters[d0] = iter - 1;
		if (objective.equals(ObjectiveFunction.MUL_KER_PCT)){	
			System.out.println("SEARCH DONE obj_val " + min_obj_val + " iters " + (iter - 1));
		}
	}

	private HashMap<Integer, int[]> setupIndicies(HashMap<Integer, int[]> legal_pairs, int[] i_Ts, int[] i_Cs) {
		HashMap<Integer, int[]> Ts_to_Cs = new HashMap<Integer, int[]>();
		if (legal_pairs == null) {
			for (int i_T : i_Ts){
				Ts_to_Cs.put(i_T, i_Cs);
			}
		} else {
			Ts_to_Cs = legal_pairs;
		}
		return Ts_to_Cs;
	}

	private void createScaledXstd() {
		Xscaled = new double[X.length][];
		int p = X[0].length;
		for (int i = 0; i < X.length; i++){
			for (int j = 0; j < p; j++){
				if (Xscaled[i] == null){
					Xscaled[i] = new double[p];
				}
				Xscaled[i][j] = X[i][j] / nT;
			}			
		}
	}

	private void updateAvgVec(double[] avg_vec, int i_remove, int i_add, int nT) {
		double[] obs_to_remove = Xscaled[i_remove];
		double[] obs_to_add = Xscaled[i_add];
		
		for (int j = 0; j < obs_to_add.length; j++){
			avg_vec[j] = avg_vec[j] - obs_to_remove[j] + obs_to_add[j];
		}		
	}
}
