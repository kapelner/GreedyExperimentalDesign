package ExperimentalDesign;

import java.util.concurrent.atomic.AtomicInteger;

public abstract class MultipleSearchExperimentalDesigns extends AllExperimentalDesigns {

	protected int max_designs;
	protected int[][] ending_indicTs;	
	protected Double[] objective_vals;	
	protected Integer[] num_iters;
	protected AtomicInteger num_completed;
	
	public void beginSearch(){
		super.beginSearch();
		
		num_completed = new AtomicInteger(0);
		num_iters = new Integer[max_designs];
		ending_indicTs = new int[max_designs][n];
		objective_vals = new Double[max_designs];
	}

	
	public void setMaxDesigns(int max_designs){
		this.max_designs = max_designs;
//		System.out.println("max_designs " + this.max_designs);
	}
	
	public int[] getNumIters(){
		int[] num_iters = new int[num_completed.get()];
		for (int i = 0; i < num_iters.length; i++){
			num_iters[i] = this.num_iters[i];
		}
		return num_iters;
	}
	
	public double[] getObjectiveVals(){
//		System.out.println("getObjectiveVals num_completed: " + num_completed);
		double[] objective_vals = new double[num_completed.get()];
		for (int i = 0; i < objective_vals.length; i++){
			objective_vals[i] = (this.objective_vals[i] == null) ? 0 : this.objective_vals[i];
		}
		return objective_vals;
	}
	
	public int[][] getEndingIndicTs(int index){ //stupid R
		int[] indicies = {index};
		return getEndingIndicTs(indicies);
	}
	
	public int[][] getEndingIndicTs(int[] indicies){
		int[][] ending_indicTs = new int[indicies.length][n];
		for (int i = 0; i < indicies.length; i++){
			ending_indicTs[i] = this.ending_indicTs[indicies[i]];
		}
		return ending_indicTs;
	}
	
	public int[][] getEndingIndicTs(){
		return ending_indicTs;
	}
	
	public int progress(){
		return num_completed.get();
	}

}
