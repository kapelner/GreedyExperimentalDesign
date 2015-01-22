package GreedyExperimentalDesign;

import no.uib.cipr.matrix.DenseVector;

public abstract class ObjectiveFunction {
	protected DenseVector XTbar;
	protected DenseVector XCbar;

	public abstract double calc();

	public void setXTbar(DenseVector XTbar){
//		System.out.println("XTbar: " + Tools.StringJoin(XTbar.getData()));
		this.XTbar = XTbar;
	}
	
	public void setXCbar(DenseVector XCbar){
//		System.out.println("XCbar: " + Tools.StringJoin(XCbar.getData()));
		this.XCbar = XCbar;
	}
}
