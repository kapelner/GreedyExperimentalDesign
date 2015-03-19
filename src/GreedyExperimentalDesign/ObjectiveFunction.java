package GreedyExperimentalDesign;

public abstract class ObjectiveFunction {
	protected double[] XTbar;
	protected double[] XCbar;

	public abstract double calc(boolean debug_mode);

	public void setXTbar(double[] XTbar){
//		System.out.println("XTbar: " + Tools.StringJoin(XTbar.getData()));
		this.XTbar = XTbar;
	}
	
	public void setXCbar(double[] XCbar){
//		System.out.println("XCbar: " + Tools.StringJoin(XCbar.getData()));
		this.XCbar = XCbar;
	}
}
