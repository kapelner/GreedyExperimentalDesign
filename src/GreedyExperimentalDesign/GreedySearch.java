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
	private boolean verbose;
	private boolean use_gpu;

	public GreedySearch(
		int nT,
		double[][] X, 
		double[][] Sinvmat, 
		ArrayList<int[]> legal_pairs,
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
		boolean semigreedy, 
		boolean diagnostics, 
		boolean verbose,
		Integer max_iters, 
		Random r, 
		AtomicBoolean search_stopped, 
		HashMap<Integer, double[][]> Kgrams, 
		HashMap<Integer, Double> max_reduction_log_obj_vals, 
		double[] kernel_weights, 
		Double maximum_gain_scaling, 
		ArrayList<double[]> kernel_obj_values,
		boolean use_gpu
	) {

		this.nT = nT;
		this.use_gpu = use_gpu;
		
		ObjectiveFunction obj_fun = null;
		this.verbose = verbose;
		int p = 0;
		int n = 0;
		if (objective.equals(ObjectiveFunction.KER)){
			obj_fun = new KernelObjective(Kgram);	
			n = Kgram.length;
			((KernelObjective)obj_fun).setW(indicT);
			((KernelObjective)obj_fun).setUseGpu(use_gpu);
		} else if (objective.equals(ObjectiveFunction.MUL_KER_PCT)) {
			obj_fun = new MultipleKernelObjectiveFunction(Kgrams, max_reduction_log_obj_vals, kernel_weights, maximum_gain_scaling, kernel_obj_values, verbose);	
			n = Kgrams.get(0).length; 
			((MultipleKernelObjectiveFunction)obj_fun).setW(indicT);
			((MultipleKernelObjectiveFunction)obj_fun).setInitialObjVals();
			((MultipleKernelObjectiveFunction)obj_fun).setUseGpu(use_gpu);
			if (diagnostics) {
				((MultipleKernelObjectiveFunction)obj_fun).resetKernelSum(); 
				((MultipleKernelObjectiveFunction)obj_fun).calcKernelObjDiagnostics();
			}
			((MultipleKernelObjectiveFunction)obj_fun).resetKernelSum(); 
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
		
		Double obj_val = null;		
		int[] i_Ts = Tools.findIndicies(indicT, nT, 1);
		int[] i_Cs = Tools.findIndicies(indicT, n - nT, 0);

		double[] avg_Ts = null;
		double[] avg_Cs = null;
		if (obj_fun instanceof SimpleAverageObjectiveFunction){
			ArrayList<double[]> XT = Tools.subsetMatrix(X, i_Ts); 
			ArrayList<double[]> XC = Tools.subsetMatrix(X, i_Cs); 
			avg_Ts = Tools.colAvg(XT, p);
			avg_Cs = Tools.colAvg(XC, p);
			((SimpleAverageObjectiveFunction)obj_fun).setXTbar(avg_Ts);
			((SimpleAverageObjectiveFunction)obj_fun).setXCbar(avg_Cs);
		}
		
		int iter = 0;
		double min_obj_val = obj_fun.calc(false);

		OptimalExperimentalDesign.WebGpuPanama.GpuGreedySearchSession gpuSession = null;
		if (use_gpu && objective.equals(ObjectiveFunction.MAHAL) && legal_pairs == null) {
			double[] XFlat = new double[n * p];
			for (int i = 0; i < n; i++) for (int j = 0; j < p; j++) XFlat[i * p + j] = X[i][j];
			double[] SinvFlat = new double[p * p];
			for (int i = 0; i < p; i++) for (int j = 0; j < p; j++) SinvFlat[i * p + j] = Sinvmat[i][j];
			gpuSession = OptimalExperimentalDesign.WebGpuPanama.createGreedySearchSession(n, p, XFlat, SinvFlat, i_Ts, i_Cs);
		}

		while (true){
			iter++;
			int[] switched_pair_1 = new int[2];
			int[] switched_pair_2 = new int[2];
			int[] indicTmin = null;
			
			if (legal_pairs == null) {
				if (gpuSession != null) {
					try {
						double[] switch_obj_vals = gpuSession.runIteration(nT, avg_Ts, avg_Cs);
						int idx_gpu = 0;
						for (int i_T : i_Ts) {
							for (int i_C : i_Cs) {
								double obj_val_gpu = switch_obj_vals[idx_gpu++];
								if (obj_val_gpu < min_obj_val) {
									min_obj_val = obj_val_gpu;
									indicTmin = indicT.clone();
									indicTmin[i_T] = 0;
									indicTmin[i_C] = 1;
									if (diagnostics) {
										switched_pair_1[0] = i_T;
										switched_pair_1[1] = i_C;
									}
								}
							}
						}
					} catch (Throwable t) {
						throw new RuntimeException(t);
					}
				} else {
					// Save kernel sums before the inner loop so each candidate
					// is evaluated from the same baseline (setSwitch mutates w in place).
					double saved_kernel_sum = 0;
					double[] saved_multi_kernel_sums = null;
					if (objective.equals(ObjectiveFunction.KER)) {
						saved_kernel_sum = ((KernelObjective)obj_fun).getRunningKernelSum();
					} else if (objective.equals(ObjectiveFunction.MUL_KER_PCT)) {
						saved_multi_kernel_sums = ((MultipleKernelObjectiveFunction)obj_fun).getRunningKernelSums();
					}
					indices_loop:
					for (int i_T : i_Ts){
						for (int i_C : i_Cs){
							if (objective.equals(ObjectiveFunction.KER)){
								((KernelObjective)obj_fun).setSwitch(i_T, i_C);
							} else if (objective.equals(ObjectiveFunction.MUL_KER_PCT)) {
								((MultipleKernelObjectiveFunction)obj_fun).setSwitch(i_T, i_C);
							} else {
								updateAvgVec(avg_Ts, i_T, i_C, nT);
								updateAvgVec(avg_Cs, i_C, i_T, n - nT);
							}
							obj_val = obj_fun.calc(false);
							boolean is_new_min = obj_val < min_obj_val;
							if (is_new_min){
								indicTmin = indicT.clone();
								indicTmin[i_T] = 0;
								indicTmin[i_C] = 1;
								min_obj_val = obj_val;
								if (diagnostics){
									switched_pair_1[0] = i_T;
									switched_pair_1[1] = i_C;
								}
							}
							// Undo switch so each candidate is evaluated from the same
							// baseline and indicT is never left with a wrong count of 1s.
							if (objective.equals(ObjectiveFunction.KER)){
								((KernelObjective)obj_fun).restoreW(i_T, i_C);
								((KernelObjective)obj_fun).setRunningKernelSum(saved_kernel_sum);
							} else if (objective.equals(ObjectiveFunction.MUL_KER_PCT)) {
								((MultipleKernelObjectiveFunction)obj_fun).restoreKernelSumsAndW(i_T, i_C, saved_multi_kernel_sums);
							} else if (obj_fun instanceof SimpleAverageObjectiveFunction){
								updateAvgVec(avg_Ts, i_C, i_T, nT);
								updateAvgVec(avg_Cs, i_T, i_C, n - nT);
							}
							if (is_new_min && semigreedy) break indices_loop;
						}
					}
				}			
			} else {
				indices_loop:
				for (int m1 = 0; m1 < (legal_pairs.size() - 1); m1++){
					for (int m2 = m1 + 1; m2 < legal_pairs.size(); m2++){
						int[] pair1 = legal_pairs.get(m1);
						int[] pair2 = legal_pairs.get(m2);
						int i_T_1 = (indicT[pair1[0]] == 1) ? pair1[0] : pair1[1];
						int i_C_1 = (indicT[pair1[0]] == 1) ? pair1[1] : pair1[0];
						int i_T_2 = (indicT[pair2[0]] == 1) ? pair2[0] : pair2[1];
						int i_C_2 = (indicT[pair2[0]] == 1) ? pair2[1] : pair2[0];
						
						if (objective.equals(ObjectiveFunction.KER)){	
							((KernelObjective)obj_fun).setSwitch(i_T_1, i_C_1);	
							((KernelObjective)obj_fun).setSwitch(i_T_2, i_C_2);						
						} else if (objective.equals(ObjectiveFunction.MUL_KER_PCT)) {
							((MultipleKernelObjectiveFunction)obj_fun).setSwitch(i_T_1, i_C_1);	
							((MultipleKernelObjectiveFunction)obj_fun).setSwitch(i_T_2, i_C_2);		
						} else {
							updateAvgVec(avg_Ts, i_T_1, i_C_1, nT);
							updateAvgVec(avg_Ts, i_T_2, i_C_2, nT);
							updateAvgVec(avg_Cs, i_C_1, i_T_1, n - nT);
							updateAvgVec(avg_Cs, i_C_2, i_T_2, n - nT);
						}
						obj_val = obj_fun.calc(false);
						if (obj_val < min_obj_val){
							indicTmin = indicT.clone();
							indicTmin[i_T_1] = 0; indicTmin[i_C_1] = 1;
							indicTmin[i_T_2] = 0; indicTmin[i_C_2] = 1;
							min_obj_val = obj_val;
							if (diagnostics){
								switched_pair_1[0] = i_T_1; switched_pair_1[1] = i_C_1;
								switched_pair_2[0] = i_T_2; switched_pair_2[1] = i_C_2;
							}
							if (semigreedy) break indices_loop;
						}
						if (obj_fun instanceof SimpleAverageObjectiveFunction){
							updateAvgVec(avg_Ts, i_C_1, i_T_1, nT);
							updateAvgVec(avg_Ts, i_C_2, i_T_2, nT);
							updateAvgVec(avg_Cs, i_T_1, i_C_1, n - nT);	
							updateAvgVec(avg_Cs, i_T_2, i_C_2, n - nT);	
						}
					}	
				}
			}

			if (diagnostics && indicTmin != null){
				switched_pairs.add(switched_pair_1);
				if (legal_pairs != null) switched_pairs.add(switched_pair_2);
				min_obj_val_by_iteration.add(min_obj_val);
			}
			
			if (indicTmin == null) break;
			indicT = indicTmin;
			i_Ts = Tools.findIndicies(indicT, nT, 1);
			i_Cs = Tools.findIndicies(indicT, n - nT, 0);

			if (objective.equals(ObjectiveFunction.KER)){	
				((KernelObjective)obj_fun).resetKernelSum();
				((KernelObjective)obj_fun).setW(indicT);
				obj_fun.calc(false);
			} else if (objective.equals(ObjectiveFunction.MUL_KER_PCT)){	
				((MultipleKernelObjectiveFunction)obj_fun).resetKernelSum();
				((MultipleKernelObjectiveFunction)obj_fun).setW(indicT);
				obj_fun.calc(false);
			}

			if (max_iters != null && max_iters == iter) break;
			if (search_stopped.get()) break;
		}	
		
		if (gpuSession != null) gpuSession.close();
		for (int i = 0; i < indicT.length; i++) ending_indicT[i] = indicT[i];
		objective_vals[d0] = min_obj_val;
		num_iters[d0] = iter - 1;
	}

	private void createScaledXstd() {
		Xscaled = new double[X.length][];
		int p = X[0].length;
		for (int i = 0; i < X.length; i++){
			Xscaled[i] = new double[p];
			for (int j = 0; j < p; j++) Xscaled[i][j] = X[i][j] / nT;
		}
	}

	private void updateAvgVec(double[] avg_vec, int i_remove, int i_add, int nT) {
		double[] obs_to_remove = Xscaled[i_remove];
		double[] obs_to_add = Xscaled[i_add];
		for (int j = 0; j < obs_to_add.length; j++) avg_vec[j] = avg_vec[j] - obs_to_remove[j] + obs_to_add[j];
	}
}
