package ObjectiveFunctions;

import java.util.ArrayList;
import java.util.HashMap;

import ExperimentalDesign.Tools;

//import org.apache.commons.math3.random.EmpiricalDistribution;

public class MultipleKernelObjectiveFunction extends ObjectiveFunction {


	
//	private HashMap<Integer, double[][]> Kgrams;
//	private HashMap<Integer, double[]> obj_val_kernel_library;
	private double[] kernel_weights;
	private int m;
	private KernelObjective[] kernel_objective_functions;
//	private EmpiricalDistribution[] kernel_ecdfs;
	private double maximum_gain_scaling;
	private HashMap<Integer, Double> max_reduction_log_obj_vals;
	private ArrayList<double[]> kernel_obj_values;

	public MultipleKernelObjectiveFunction(
			HashMap<Integer, double[][]> Kgrams,	
//			HashMap<Integer, double[]> obj_val_kernel_library, 
			HashMap<Integer, Double> max_reduction_log_obj_vals, 
			double[] kernel_weights, 
			double maximum_gain_scaling, 
			ArrayList<double[]> kernel_obj_values
		) {
//		this.Kgrams = Kgrams;
//		this.obj_val_kernel_library = obj_val_kernel_library;
		this.kernel_weights = kernel_weights;
		this.max_reduction_log_obj_vals = max_reduction_log_obj_vals;
		this.maximum_gain_scaling = maximum_gain_scaling;
		this.kernel_obj_values = kernel_obj_values;
		m = Kgrams.size();
		
		//now we create sub-objective functions
		kernel_objective_functions = new KernelObjective[m];
		for (int i_k = 0; i_k < m; i_k++) {
			kernel_objective_functions[i_k] = new KernelObjective(Kgrams.get(i_k));
		}
		System.out.println("MultipleKernelObjectiveFunction init m = " + m + " maximum_gain_scaling =  " + maximum_gain_scaling + " kernel_weights = " + Tools.StringJoin(kernel_weights));
		//now cache the empirical distribution object
//		for (int i_k = 0; i_k < m; i_k++) {
//			kernel_ecdfs[i_k] = new EmpiricalDistribution();
//			kernel_ecdfs[i_k].load(obj_val_kernel_library.get(i_k));
//		}
	}

	@Override
	public double calc(boolean debug_mode) {
		double obj_val = 0;
		double[] obj_vals = new double[m];
		for (int i_k = 0; i_k < m; i_k++) {
//			System.out.print("    i_k " + (i_k + 1));	
			double obj_val_k = (kernel_weights[i_k] * pct_off_from_best(i_k));
			obj_vals[i_k] = obj_val_k;
			obj_val += obj_val_k;
		}
//		System.out.println("    aggregate objval: " + obj_val + "\n");
		kernel_obj_values.add(obj_vals);
		return obj_val;
	}

	private double pct_off_from_best(int i_k) {
//		System.out.print(" max_reduction " + String.format("%.4g", max_reduction_log_obj_vals.get(i_k)));	
//		kernel_ecdfs[i_k].cumulativeProbability(x);
		double ret = 1 - kernel_objective_functions[i_k].log10_i_over_current_obj_val() / (max_reduction_log_obj_vals.get(i_k) * maximum_gain_scaling);
//		System.out.print(" pct_off_from_best " + String.format("%.4g", ret) + "\n");
		
		return ret; 
	}

	public void setW(int[] indicT) {
		for (int i_k = 0; i_k < m; i_k++) {
			kernel_objective_functions[i_k].setW(indicT);
		}
	}

	public void setSwitch(int t, int c) {
//		System.out.println("    switch " + t + " <> " + c);
		for (int i_k = 0; i_k < m; i_k++) {
			kernel_objective_functions[i_k].setSwitch(t, c);
		}
	}

	public void resetKernelSum() {
		for (int i_k = 0; i_k < m; i_k++) {
			kernel_objective_functions[i_k].resetKernelSum();
		}
	}

	public void setInitialObjVals() {
		for (int i_k = 0; i_k < m; i_k++) {
			kernel_objective_functions[i_k].setInitialObjVal();
//			System.out.println("  i_k " + i_k + " initial objval: " + kernel_objective_functions[i_k].running_kernel_sum);	
		}
	}

}
