package GreedyExperimentalDesign;

import no.uib.cipr.matrix.DenseVector;

public abstract class ObjectiveFunction {
	protected DenseVector XTbar;
	protected DenseVector XCbar;

	public abstract double calc();

	public void setXTbar(DenseVector XTbar){
		this.XTbar = XTbar;
	}
	
	public void setXCbar(DenseVector XCbar){
		this.XCbar = XCbar;
	}
}
