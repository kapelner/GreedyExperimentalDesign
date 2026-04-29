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

	// Read-only O(n) proposal evaluation — non-mutating equivalent of pre-GPU calc().
	// obj_new = obj_old + 4 * sum_{l != i_T, i_C} v_l * (K[i_C][l] - K[i_T][l])
	public double calcProposal(int i_T, int i_C) {
		double delta = 0.0;
		for (int l = 0; l < n; l++) {
			if (l == i_T || l == i_C) continue;
			delta += (2 * w[l] - 1) * (Kgram[i_C][l] - Kgram[i_T][l]);
		}
		return running_kernel_sum + 4 * delta;
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
