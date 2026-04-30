package ObjectiveFunctions;

public class KernelObjective extends ObjectiveFunction {

	private double[][] Kgram;
	private int[] w;
	private int n;
	public double running_kernel_sum;
	public double initial_obj_val;
	
	public KernelObjective(double[][] Kgram) {
		this.Kgram = Kgram;
		this.n = Kgram.length;
	}

	@Override
	public double calc(boolean b) {
		return running_kernel_sum;
	}

	public void setW(int[] w) {
		this.w = w;
		running_kernel_sum = 0;
		calcCpu();
	}

	// Read-only O(n) proposal evaluation — non-mutating equivalent of setSwitch() followed by calc().
	// v_new = v + delta_v where delta_v_iT = -2 and delta_v_iC = 2
	// obj_new = v_new^T K v_new = (v + delta_v)^T K (v + delta_v) = v^T K v + 2 v^T K delta_v + delta_v^T K delta_v
	public double calcProposal(int i_T, int i_C) {
		// delta_v^T K delta_v = (-2)^2 * K_iT,iT + (2)^2 * K_iC,iC + 2 * (-2) * (2) * K_iT,iC
		// = 4 * K_iT,iT + 4 * K_iC,iC - 8 * K_iT,iC
		double term3 = 4 * Kgram[i_T][i_T] + 4 * Kgram[i_C][i_C] - 8 * Kgram[i_T][i_C];
		
		// 2 * v^T K delta_v = 2 * sum_j v_j * (K delta_v)_j = 2 * sum_j v_j * (K_j,iT * (-2) + K_j,iC * (2))
		// = 4 * sum_j v_j * (K_j,iC - K_j,iT)
		double term2 = 0;
		for (int j = 0; j < n; j++) {
			double v_j = 2 * w[j] - 1;
			term2 += v_j * (Kgram[i_C][j] - Kgram[i_T][j]);
		}
		return running_kernel_sum + 4 * term2 + term3;
	}
    
    private void calcCpu() {
        running_kernel_sum = 0;
        for (int i = 0; i < n; i++) {
            double vi = 2 * w[i] - 1;
            for (int j = 0; j < n; j++) {
                double vj = 2 * w[j] - 1;
                running_kernel_sum += vi * Kgram[i][j] * vj;
            }
        }
    }

	public void setSwitch(int i_T, int i_C) {
		//we are removing treatment from i_T and adding it to i_C
		//this means vi goes from 1 to -1 and vj goes from -1 to 1
		//delta = (v_new - v_old)_i * K_i,i * (v_new - v_old)_i + 2 * (v_new - v_old)_i * sum_{j != i} K_i,j * v_j
		
		//update for i_T
		double v_old_iT = 1;
		double v_new_iT = -1;
		double diff_iT = v_new_iT - v_old_iT; // -2
		
		double delta = diff_iT * Kgram[i_T][i_T] * diff_iT;
		for (int j = 0; j < n; j++) {
			if (j == i_T) {
				continue;
			}
			double v_j = 2 * w[j] - 1;
			delta += 2 * diff_iT * Kgram[i_T][j] * v_j;
		}
		running_kernel_sum += delta;
		w[i_T] = 0;
		
		//update for i_C
		double v_old_iC = -1;
		double v_new_iC = 1;
		double diff_iC = v_new_iC - v_old_iC; // 2
		
		delta = diff_iC * Kgram[i_C][i_C] * diff_iC;
		for (int j = 0; j < n; j++) {
			if (j == i_C) {
				continue;
			}
			double v_j = 2 * w[j] - 1;
			delta += 2 * diff_iC * Kgram[i_C][j] * v_j;
		}
		running_kernel_sum += delta;
		w[i_C] = 1;
	}

	public void resetKernelSum() {
		setW(w);
	}

	public void setInitialObjVal() {
		initial_obj_val = running_kernel_sum;
	}

	public double log10_i_over_current_obj_val() {
		return Math.log10(initial_obj_val / running_kernel_sum);
	}

	public double getRunningKernelSum() {
		return running_kernel_sum;
	}

	public void setRunningKernelSum(double sum) {
		running_kernel_sum = sum;
	}

	public void restoreW(int i_T, int i_C) {
		w[i_T] = 1;
		w[i_C] = 0;
	}

}
