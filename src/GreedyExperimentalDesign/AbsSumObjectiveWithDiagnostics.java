package GreedyExperimentalDesign;

import java.util.ArrayList;

public class AbsSumObjectiveWithDiagnostics extends ObjectiveFunction {
	
	private ArrayList<double[]> xbardiffjs_by_iteration;
	
	public AbsSumObjectiveWithDiagnostics(ArrayList<double[]> xbardiffjs_by_iteration){
		this.xbardiffjs_by_iteration = xbardiffjs_by_iteration;
	}

	public ArrayList<double[]> getXbardiffjs_by_iteration() {
		return xbardiffjs_by_iteration;
	}

	@Override
	public double calc(boolean debug_mode) {
		if (debug_mode){
			System.out.println("XTbar: " + Tools.StringJoin(XTbar, ","));
			System.out.println("XCbar: " + Tools.StringJoin(XCbar, ","));			
		}
		int p = XTbar.length;
		double[] xbardiffjs = new double[p];
		double abs_sum = 0;
		for (int j = 0; j < p; j++){
			double diff = XTbar[j] - XCbar[j];
			xbardiffjs[j] = diff;
			abs_sum += (diff < 0.0 ? -diff : diff); //faster than Math.abs accd to JProfiler
		}		
		xbardiffjs_by_iteration.add(xbardiffjs);
		return abs_sum;
	}

}
