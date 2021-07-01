package ObjectiveFunctions;

import java.util.HashMap;

public class KernelObjective extends ObjectiveFunction {


	private double[][] Kgram;
//	private HashMap<Integer, HashMap<Integer, Double>> qcds;
//	private double[][] qcds;
	protected int[] w;
	private int n;
	protected Double running_kernel_sum;
	private Integer t; //the index of the new treatment
	private Integer c; //the index of the new control

	public KernelObjective(double[][] Kgram) {
		this.Kgram = Kgram;
		this.n = Kgram.length;
//		qcds = new HashMap<Integer, HashMap<Integer, Double>>();
//		qcds = new double[n][n];
	}

	@Override
	public double calc(boolean debug_mode) {
		//we've started with a new vector so need to do the full calculation once
		if (running_kernel_sum == null) {
//			System.out.println("running_kernel_sum == null");
			fullQuadraticFormCalculationAndCache(debug_mode);
			if (t == null) { //for diagnostics only
				return running_kernel_sum;
			}
		}
		double qcd = 0.0;
//		qcds[t][c] = 0.0;
		for (int l = 0; l < n; l++) {
			if (l == t || l == c) {
				continue;
			}
			qcd += (w[l] * (Kgram[t][l] - Kgram[c][l]));	
//			qcds[t][c] += (w[l] * (Kgram[t][l] - Kgram[c][l]));	
		}
		
		//cache
//		if (qcds.get(t) == null) {
//			qcds.put(t, new HashMap<Integer, Double>());
//		}
//		qcds.get(t).put(c, qcd);
		//return
		return running_kernel_sum - 4 * qcd;
//		return running_kernel_sum - 4 * qcds[t][c];
	}
	
//	public void setPermanentSwitch(int t, int c) {
//		//set the running kernel sum
////		setSwitch(t, c);
////		running_kernel_sum = calc(false);
////		running_kernel_sum -= 4 * qcds.get(t).get(c);
//		running_kernel_sum -= 4 * qcds[t][c];
//		//finally, reset cache and temp switch values
//		this.t = null;
//		this.c = null;
////		qcds = new HashMap<Integer, HashMap<Integer, Double>>();
//		qcds = new double[n][n];
//	}
	
	public void resetKernelSum() {
		running_kernel_sum = null;
	}
	
	private void fullQuadraticFormCalculationAndCache(boolean debug_mode) {
		this.t = null;
		this.c = null;
//		qcds = new HashMap<Integer, HashMap<Integer, Double>>();
//		qcds = new double[n][n];
//		System.out.println("fullQuadraticFormCalculationAndCache Kgram = " + Kgram + " w = " + w + " running_kernel_sum = " + running_kernel_sum);
		running_kernel_sum = 0.0;
		
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				running_kernel_sum += Kgram[i][j] * w[i] * w[j];
			}
		}
		//we're leaving out this constant
		//kernel_sum *= 4 / n^2; //this is the constant in Eq 4.2 of Kallus (2018)
	
		if (debug_mode){
			System.out.println("kernel_sum: " + running_kernel_sum);			
		}
	}
	
	public void setSwitch(int t, int c) {
		this.t = t;
		this.c = c;
	}

	//binary => 1/-1. This setting is done once
	public void setW(int[] indicT) {
		int n = indicT.length;
		this.w = new int[n];
		for (int i = 0; i < n; i++) {
			this.w[i] = (indicT[i] == 1 ? 1 : -1);
		}
	}

}
