package GreedyExperimentalDesign;

import no.uib.cipr.matrix.DenseMatrix;

public class PropMahalObjective extends ObjectiveFunction {

	private DenseMatrix Sinvmat;

	public PropMahalObjective(DenseMatrix Sinvmat) {
		this.Sinvmat = Sinvmat;
	}

	@Override
	public double calc() {
		//as.numeric(t(X_T_bar_minus_X_C_bar) %*% inv_cov_X %*% X_T_bar_minus_X_C_bar)
		int p = XTbar.size();
		DenseMatrix X_T_bar_minus_X_C_bar = new DenseMatrix(p, 1);
		for (int j = 0; j < p; j++){
			X_T_bar_minus_X_C_bar.set(j, 0, XTbar.get(j) - XCbar.get(j));
		}
		DenseMatrix temp = new DenseMatrix(1, p);
		(X_T_bar_minus_X_C_bar.transpose()).mult(Sinvmat, temp);
		DenseMatrix temp2 = new DenseMatrix(1, 1);
		temp.mult(X_T_bar_minus_X_C_bar, temp2);
		return temp2.get(0, 0); //it's a scalar at the end
	}

}
