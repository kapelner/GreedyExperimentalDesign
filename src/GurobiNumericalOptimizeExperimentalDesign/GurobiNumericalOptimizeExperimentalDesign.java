package GurobiNumericalOptimizeExperimentalDesign;

import ExperimentalDesign.AllExperimentalDesigns;
import ObjectiveFunctions.ObjectiveFunction;

import org.ejml.simple.SimpleMatrix;
import gurobi.*;

public class GurobiNumericalOptimizeExperimentalDesign extends AllExperimentalDesigns {

	/** the value to be returned after optimization */	
	private int[] indicator_T;
	/** how long can the optimizer take? */	
	private Double time_limit_min;
	/** how many nodes can the optimizer explore? */	
	private Integer node_limit;
	
	//running the Java as standalone is for debug purposes ONLY!!!
	public static void main(String[] args) {	

		GurobiNumericalOptimizeExperimentalDesign gnuoed = new GurobiNumericalOptimizeExperimentalDesign();
		//set seed here for reproducibility during debugging
		gnuoed.rand_obj.setSeed(1984);

		int n = 100;
		int p = 10;
		try {
			gnuoed.setNandP(n, p);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		for (int i = 0; i < n; i++){
			double[] x_i = new double[p];
			for (int j = 0; j < p; j++){
				x_i[j] = gnuoed.rand_obj.nextDouble();
			}
			gnuoed.setDataRow(i, x_i);
		}
		gnuoed.setNumCores(3);
		gnuoed.setWait();
		gnuoed.beginSearch();
	}


	
	public void beginSearch() {
		super.beginSearch();
		
		//bracha: the data "X" is the data matrix and "Sinv" is the inverse sample var-cov matrix"
				
	    GRBEnv env = null;    	
	    	
	        
			try {
				env = new GRBEnv("gurobi_numerical_optimization_via_R_package_Gree.log");
			} catch (GRBException e) {
				System.err.println("Gurobi error when creating the environment. Error code: " + e.getErrorCode());
				e.printStackTrace();
			}
	        try {
				env.set(GRB.IntParam.Threads, num_cores);
			} catch (GRBException e) {
				System.err.println("Gurobi error when setting the number of cores. Error code: " + e.getErrorCode());
				e.printStackTrace();
			}
	        
	        GRBModel model = null;
			try {
				model = new GRBModel(env);
			} catch (GRBException e) {
				System.err.println("Gurobi error when creating the model. Error code: " + e.getErrorCode());
				e.printStackTrace();
			}
			if (time_limit_min != null){
				 try {
					model.set(GRB.DoubleParam.TimeLimit, time_limit_min * 60);
				} catch (GRBException e) {
					System.err.println("Gurobi error when setting the time limit. Error code: " + e.getErrorCode());
					e.printStackTrace();
				}
			}
			if (node_limit != null){
				 try {
					model.set(GRB.DoubleParam.NodeLimit, node_limit);
				} catch (GRBException e) {
					System.err.println("Gurobi error when setting the time limit. Error code: " + e.getErrorCode());
					e.printStackTrace();
				}
			}
	       
	       
	            	
	    	SimpleMatrix Xsm = new SimpleMatrix(X);
	        SimpleMatrix Sinvsm = new SimpleMatrix(Sinv);
	    	SimpleMatrix XSinvXt = Xsm.mult(Sinvsm).mult(Xsm.transpose());
//	    	System.out.println("s: "+ Sinvsm.get(0,0));
//	    	Xsm = null;
	    	
//	    	System.out.println("sinv: " + Sinvsm.get(0,0) + " n: "+n+" p : "+p);
	    	

	        // Create variables
	    	GRBVar[] indicator_T_gur = new GRBVar[n];
	    	for (int i = 0; i < n; i++) {
	    		try {
    	    		indicator_T_gur[i] = model.addVar(0.0, 1.0, 0.0, GRB.BINARY, Integer.toString(i));
    	    	}
				catch (GRBException e) {
					System.err.println("Gurobi error when setting the time limit for observation #" + i + " Error code: " + e.getErrorCode());
					e.printStackTrace();
				}
	    	}
	    	
	    	//Setting objective matrix
	    	GRBQuadExpr obj = new GRBQuadExpr();
	    	for (int i = 0; i < n; i++) {
	    		for (int j = 0; j < n; j++) {
	    			obj.addTerm(1000000 * XSinvXt.get(i, j), indicator_T_gur[i], indicator_T_gur[j]);	    			
	    		}
	    	}
	    	try {
				model.setObjective(obj);
			} catch (GRBException e) {
				System.err.println("Gurobi error when setting the objective function in the model. Error code: " + e.getErrorCode());
				e.printStackTrace();
			}


	        // Add constraint: sum of vars equal to n/2 (equal num 1's and 0's)
	    	GRBLinExpr expr = new GRBLinExpr();
	    	for (int i = 0; i < n; i++) {
	    		 expr.addTerm(1.0, indicator_T_gur[i]);
	    	}
	    	
	        try {
				model.addConstr(expr, GRB.EQUAL, n / 2, "c0");
			} catch (GRBException e) {
				System.err.println("Gurobi error when setting the constraint of equal treatments and equal controls. Error code: " + e.getErrorCode());
				e.printStackTrace();
			} //ensures equal number of treatments and controls

	        // Optimize model
	        try {
				model.optimize();
			} catch (GRBException e) {
				System.err.println("Gurobi error when running the optimization algorithm. Error code: " + e.getErrorCode());
				e.printStackTrace();
			}
	        
//	        for (int i = 0; i < n; i++) {
//	        	System.out.println(indicator_T_gur[i].get(GRB.StringAttr.VarName)
//                        + " " +indicator_T_gur[i].get(GRB.DoubleAttr.X));
//	    	}

//	        System.out.println("Obj: " + model.get(GRB.DoubleAttr.ObjVal));
	        

	        
	        
	        //convert Gurobi indicator to a int vector
	    	indicator_T = new int[n];
	    	for (int i = 0; i < n; i++) {
	    		indicator_T[i] = -99;
	    	}
	    	for (int i = 0; i < n; i++) {
	    		try {
					indicator_T[i] = (int)indicator_T_gur[i].get(GRB.DoubleAttr.X);
	    			System.out.println(indicator_T_gur[i].get(GRB.StringAttr.VarName)
	                         + " " +indicator_T_gur[i].get(GRB.DoubleAttr.X));
				} catch (GRBException e) {
					System.err.println("Gurobi error when extracting the solution for vector element #" + i + " Error code: " + e.getErrorCode());
					e.printStackTrace();
				}
	    	}

	    	// Dispose of model and environment
	        model.dispose();
	        try {
				env.dispose();
			} catch (GRBException e) {
				System.err.println("Gurobi error when disposing of the environment. Error code: " + e.getErrorCode());
				e.printStackTrace();
			}
	        
//	      } catch (GRBException e) {
//	        System.out.println("Error code: " + e.getErrorCode() + ". " +
//	                           e.getMessage());
//	        
//	      }
	}

	public void setTimeLimitMin(double time_limit_min) {
		this.time_limit_min = time_limit_min;
	}
	
	public void setTimeLimitMin(int node_limit) {
		this.node_limit = node_limit;
	}
	
	public int[] getBestIndicT() {
		return indicator_T;
	}
	
	public void setObjective(String objective) throws Exception{
		if (objective.equals(ObjectiveFunction.MAHAL)){
			throw new Exception("For Numerical Optimization, Mahalanobis is the only implemented objective function.");
		}
		this.objective = objective;
	}
}
