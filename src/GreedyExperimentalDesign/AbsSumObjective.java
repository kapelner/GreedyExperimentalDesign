package GreedyExperimentalDesign;

public class AbsSumObjective extends ObjectiveFunction {

	@Override
	public double calc() {
		double abs_sum = 0;
		for (int j = 0; j < XTbar.size(); j++){
			abs_sum += Math.abs(XTbar.get(j) - XCbar.get(j));
		}		
		return abs_sum;
	}

}
