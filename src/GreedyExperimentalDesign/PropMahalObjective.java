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
//		System.out.println("p = " + p);
		DenseMatrix X_T_bar_minus_X_C_bar = new DenseMatrix(p, 1);
		DenseMatrix X_T_bar_minus_X_C_bar_t = new DenseMatrix(1, p);
//		System.out.println("X_T_bar_minus_X_C_bar.toString() = " + X_T_bar_minus_X_C_bar.toString());
		for (int j = 0; j < p; j++){
			X_T_bar_minus_X_C_bar.set(j, 0, XTbar.get(j) - XCbar.get(j));
			X_T_bar_minus_X_C_bar_t.set(0, j, XTbar.get(j) - XCbar.get(j));
		}
		
//		System.out.println("X_T_bar_minus_X_C_bar.toString() = " + X_T_bar_minus_X_C_bar.toString());
		DenseMatrix temp = new DenseMatrix(1, p);
//		System.out.println("temp.toString() = " + temp.toString());
		X_T_bar_minus_X_C_bar_t.mult(Sinvmat, temp);
//		System.out.println("temp.toString() = " + temp.toString());
		DenseMatrix temp2 = new DenseMatrix(1, 1);
//		System.out.println("temp2.toString() = " + temp2.toString());
		temp.mult(X_T_bar_minus_X_C_bar, temp2);
//		System.out.println("temp2.toString() = " + temp2.toString());
		return temp2.get(0, 0); //it's a scalar at the end
	}

}
