package ObjectiveFunctions;

public abstract class SimpleAverageObjectiveFunction extends ObjectiveFunction {

	protected double[] XTbar;
	protected double[] XCbar;


	public void setXTbar(double[] XTbar){
//		System.out.println("XTbar: " + Tools.StringJoin(XTbar.getData()));
		this.XTbar = XTbar;
	}
	
	public void setXCbar(double[] XCbar){
//		System.out.println("XCbar: " + Tools.StringJoin(XCbar.getData()));
		this.XCbar = XCbar;
	}

}
