package ObjectiveFunctions;

import java.util.ArrayList;

public abstract class ObjectiveFunction {

	//valid objective functions
	public static final String MAHAL = "mahal_dist";
	public static final String ABS = "abs_sum_diff";
	public static final String KER = "kernel";
	public static final String MUL_KER_PCT = "added_pct_reduction";
	private static final ArrayList<String> VALID_OBJ_FUNCTIONS = new ArrayList<String>();
	static {
		VALID_OBJ_FUNCTIONS.add(MAHAL);
		VALID_OBJ_FUNCTIONS.add(ABS);
		VALID_OBJ_FUNCTIONS.add(KER);
		VALID_OBJ_FUNCTIONS.add(MUL_KER_PCT);
	};
	

	public abstract double calc(boolean debug_mode);
	
	public static boolean isValidObjFunction(String objective){
		return VALID_OBJ_FUNCTIONS.contains(objective);
	}
}
