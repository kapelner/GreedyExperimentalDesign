package ExperimentalDesign;

import org.apache.commons.math3.util.CombinatoricsUtils;

public class RandomizationMetrics {
	//all the vectors
	private int n;
	private int num_vecs;
	private int[][] ending_indicTs;
	private double[][] phats;
	private double entropy_metric;
	private double se_metric;
	
	public RandomizationMetrics(){}
	
	public void setNandR(int n, int r){
		this.n = n;
		num_vecs = r;
		ending_indicTs = new int[n][r];
	}
	
	public void setDesign(int j0, int[] indicT){
//		System.out.println("setDesign " + i0 + "  " + indicT);

		for (int i = 0; i < n; i++){
			ending_indicTs[i][j0] = indicT[i];
		}
	}
	
	public void compute(){
		phats = new double[n][n];
		//for each pair we estimate the probability the randomization
		//produces a different assignment		
		for (int i1 = 0; i1 < n - 1; i1++){
			for (int i2 = i1 + 1; i2 < n; i2++){
				int num_diff = 0;
				for (int r = 0; r < num_vecs; r++){
					num_diff += ((ending_indicTs[i1] == ending_indicTs[i2]) ? 1 : 0);
				}
				phats[i1][i2] = num_diff / num_vecs;
			}
		}
		
		//this is the probability that a random assignment is the same as another one
		double s_n = (n - 2) / ((double)(2 * n - 2));
		//number of pairs
		long num_pairs = CombinatoricsUtils.binomialCoefficient(n, 2);
		
		//calculate functions of each p_ij
		double sum_entropies = 0;
		double sum_sqd_dev = 0;
		for (int i1 = 0; i1 < n - 1; i1++){
			for (int i2 = i1 + 1; i2 < n; i2++){
				double p_hat = phats[i1][i2];
				sum_entropies += (probTimesLogProb(p_hat) + probTimesLogProb(1 - p_hat));
				sum_sqd_dev += Math.pow(p_hat - s_n, 2);
			}
		}
		
		//now calculate entropy
		double entropy_norm_factor = s_n * Math.log(s_n) + (1 - s_n) * Math.log(s_n);
		entropy_metric = 1 / ((double) num_pairs) * sum_entropies / entropy_norm_factor;
		
		//now calculate se
		double const_factor = 2 / ((double) n) * Math.sqrt((2 * n - 2) / (double)(n - 2));
		se_metric = const_factor * Math.sqrt(sum_sqd_dev);
	}

	private double probTimesLogProb(double p){
		if (p == 0){ //well known limit convention
			return 0;
		}
		return p * Math.log(p);
	}
	
	public double getRandEntropyMetric() {
		return entropy_metric;
	}

	public double getRandStdErrMetric() {
		return se_metric;
	}

	public double[][] getPhats() {
		return phats;
	}
}
