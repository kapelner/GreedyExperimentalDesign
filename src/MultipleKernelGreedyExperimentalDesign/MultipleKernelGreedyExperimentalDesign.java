package MultipleKernelGreedyExperimentalDesign;

import java.util.ArrayList;
import java.util.HashMap;

import GreedyExperimentalDesign.*;

public class MultipleKernelGreedyExperimentalDesign extends GreedyExperimentalDesign {

//	private HashMap<Integer, double[]> obj_val_kernel_library;
	private HashMap<Integer, double[][]> Kgrams;
	private double[] kernel_weights;
	private Double maximum_gain_scaling;
	private HashMap<Integer, Double> max_reduction_log_obj_vals;
	private ArrayList<ArrayList<double[]>> kernel_obj_values_by_design_num;
	
	public MultipleKernelGreedyExperimentalDesign() {
		super();
//		obj_val_kernel_library = new HashMap<Integer, double[]>();
		max_reduction_log_obj_vals = new HashMap<Integer, Double>();
		Kgrams = new HashMap<Integer, double[][]>();
	}
	
	
	public void beginSearch(){
		//initialize all data
		kernel_obj_values_by_design_num = new ArrayList<ArrayList<double[]>>(max_designs);
		for (int d = 0; d < max_designs; d++){
			kernel_obj_values_by_design_num.add(new ArrayList<double[]>());
		}	
		super.beginSearch();
	}
	
	protected GreedySearch generateIndividualSearch(int d0) {
		return new GreedySearch(
			nT,
			X,
			Sinv, 
			legal_pairs,
			null,
			starting_indicTs[d0], 
			ending_indicTs[d0], 
			switched_pairs.get(d0),
			min_obj_val_by_iterations.get(d0),
			null,
			objective_vals, 
			num_iters, 
			objective, 
			d0, 
			semigreedy, 
			diagnostics, 
			max_iters, 
			rand_obj,
			search_stopped,
			//special multiple kernel params
			Kgrams,
			max_reduction_log_obj_vals,
			kernel_weights,
			maximum_gain_scaling,
			kernel_obj_values_by_design_num.get(d0)
		);		
	}
	
	public void setSpecificKgramByRow(int k, int i0, double[] kgram_i){
		if (Kgrams.get(k) == null){
			Kgrams.put(k, new double[n][n]);
		}
		for (int j = 0; j < n; j++){
			Kgrams.get(k)[i0][j] = kgram_i[j];
		}
	}
	
//	public void setObjValLibrary(int k, double[] vals) {
//		obj_val_kernel_library.put(k, vals);
//	}
	
	public void setMaxReductionLogObjVal(int k, double val) {
		max_reduction_log_obj_vals.put(k, val);
	}
	
	public void setKernelWeights(double[] kernel_weights) {
		this.kernel_weights = kernel_weights;
	}
	
	public void setMaximumGainScaling(double maximum_gain_scaling) {
		this.maximum_gain_scaling = maximum_gain_scaling;
	}
	
	public double[][][] getObjValuesByKernel(int[] indicies){		
		double[][][] kernel_obj_vals_by_design_num = new double[indicies.length][][];
		for (int i = 0; i < indicies.length; i++){
			ArrayList<double[]> kernel_obj_values_all_iters = kernel_obj_values_by_design_num.get(indicies[i]);
			int num_iters = kernel_obj_values_all_iters.size();
			double[][] kernel_obj_values = new double[num_iters][];
			for (int j = 0; j < num_iters; j++){
				kernel_obj_values[j] = kernel_obj_values_all_iters.get(j);
			}
			kernel_obj_vals_by_design_num[i] = kernel_obj_values;
		}
		return kernel_obj_vals_by_design_num;
	}
	
}
