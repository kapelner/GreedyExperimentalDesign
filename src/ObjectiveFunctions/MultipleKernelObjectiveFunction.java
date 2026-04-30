package ObjectiveFunctions;

import java.util.ArrayList;
import java.util.HashMap;

public class MultipleKernelObjectiveFunction extends ObjectiveFunction {

	private HashMap<Integer, double[][]> Kgrams;
	private double[] kernel_weights;
	private Double maximum_gain_scaling;
	private HashMap<Integer, Double> max_reduction_log_obj_vals;
	private ArrayList<double[]> kernel_obj_values;
	private int m;
	private KernelObjective[] kernel_objective_functions;
	private boolean verbose;

	public MultipleKernelObjectiveFunction(HashMap<Integer, double[][]> Kgrams, 
			HashMap<Integer, Double> max_reduction_log_obj_vals, 
			double[] kernel_weights, Double maximum_gain_scaling, 
			ArrayList<double[]> kernel_obj_values, boolean verbose) {
		this.Kgrams = Kgrams;
		this.max_reduction_log_obj_vals = max_reduction_log_obj_vals;
		this.kernel_weights = kernel_weights;
		this.maximum_gain_scaling = maximum_gain_scaling;
		this.kernel_obj_values = kernel_obj_values;
		this.verbose = verbose;
		this.m = kernel_weights.length;
		this.kernel_objective_functions = new KernelObjective[m];
		for (int k = 0; k < m; k++) {
			this.kernel_objective_functions[k] = new KernelObjective(Kgrams.get(k));
		}
	}

	@Override
	public double calc(boolean debug_mode) {
		double obj = 0;
		for (int k = 0; k < m; k++) {
			double val = kernel_objective_functions[k].calc(debug_mode);
			if (maximum_gain_scaling != null) {
				val = 1.0 - kernel_objective_functions[k].log10_i_over_current_obj_val() / (max_reduction_log_obj_vals.get(k) * maximum_gain_scaling);
			}
			obj += kernel_weights[k] * val;
		}
		return obj;
	}

	public void setW(int[] w) {
		for (int k = 0; k < m; k++) {
			kernel_objective_functions[k].setW(w);
		}
	}

	public void setSwitch(int i_T, int i_C) {
		for (int k = 0; k < m; k++) {
			kernel_objective_functions[k].setSwitch(i_T, i_C);
		}
	}

	public void setUseGpu(boolean use_gpu) {
		this.use_gpu = use_gpu;
		for (int k = 0; k < m; k++) {
			kernel_objective_functions[k].setUseGpu(use_gpu);
		}
	}

	public void resetKernelSum() {
		for (int k = 0; k < m; k++) {
			kernel_objective_functions[k].resetKernelSum();
		}
	}

	public void calcKernelObjDiagnostics() {
		double[] res = new double[m];
		for (int k = 0; k < m; k++) {
			res[k] = kernel_objective_functions[k].running_kernel_sum;
		}
		kernel_obj_values.add(res);
	}

	public void setInitialObjVals() {
		for (int k = 0; k < m; k++) {
			kernel_objective_functions[k].setInitialObjVal();
		}
	}

	public double[] getRunningKernelSums() {
		double[] sums = new double[m];
		for (int k = 0; k < m; k++) {
			sums[k] = kernel_objective_functions[k].getRunningKernelSum();
		}
		return sums;
	}

	public double calcProposal(int i_T, int i_C) {
		double obj = 0;
		for (int k = 0; k < m; k++) {
			double val = kernel_objective_functions[k].calcProposal(i_T, i_C);
			if (maximum_gain_scaling != null) {
				val = 1.0 - Math.log10(kernel_objective_functions[k].initial_obj_val / val) / (max_reduction_log_obj_vals.get(k) * maximum_gain_scaling);
			}
			obj += kernel_weights[k] * val;
		}
		return obj;
	}


	public void restoreKernelSumsAndW(int i_T, int i_C, double[] saved_sums) {
		for (int k = 0; k < m; k++) {
			kernel_objective_functions[k].restoreW(i_T, i_C);
			kernel_objective_functions[k].setRunningKernelSum(saved_sums[k]);
		}
	}

}
