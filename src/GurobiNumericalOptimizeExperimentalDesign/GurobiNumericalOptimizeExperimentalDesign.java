package GurobiNumericalOptimizeExperimentalDesign;

import ExperimentalDesign.AllExperimentalDesigns;
import ObjectiveFunctions.ObjectiveFunction;

import org.ejml.simple.SimpleMatrix;
import gurobi.*;

public class GurobiNumericalOptimizeExperimentalDesign extends AllExperimentalDesigns {

	/** the value to be returned after optimization */	
	private int[][] indicator_T;
	/** how long can the optimizer take? */	
	private Double time_limit_min;
	/** how many nodes can the optimizer explore? */	
	private Integer node_limit;
	/** What tolerance to use? */	
	private Double tol;
	/** max number of solutions Gurobi should retain */
	private Integer max_solutions;
	/** should we turn the optimizer's screen log off? */	
	private boolean gurobi_log_off;
	/** log filename? */
	private String log_file;
	/** the Gurobi model object that will do the optimization */
	private GRBModel model;
	/** Should many solutions be close to optimal? */
	private boolean all_optimal;
	// Gurobi environment object
	private GRBEnv env;
	
	//running the Java as standalone is for debug purposes ONLY!!!
	public static void main(String[] args) {	

		GurobiNumericalOptimizeExperimentalDesign gnuoed = new GurobiNumericalOptimizeExperimentalDesign();
		//set seed here for reproducibility during debugging
		gnuoed.rand_obj.setSeed(1984);

		int n = 100;
		int p = 10;
		try {
			gnuoed.setN(n);
			gnuoed.setP(p);
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

	private void initEnv() {
		String err_msg = "Gurobi error when ";

		try {
			env = new GRBEnv();
			err_msg = "creating the environment. ";

			env.set(GRB.StringParam.LogFile, log_file);
			err_msg = "setting the log file. ";

			env.set(GRB.IntParam.Threads, num_cores);
			err_msg = "setting the number of cores. ";

			env.set(GRB.DoubleParam.MIPGapAbs, 1.00e-20);
			env.set(GRB.DoubleParam.MIPGap, 1.00e-20);
			err_msg = "setting MIP Gap. ";
			if (gurobi_log_off) {
				env.set(GRB.IntParam.LogToConsole, 0);
				err_msg = "turning off console log. ";
			}
		} catch (GRBException e) {
			System.err.println(err_msg + "Error code: " + e.getErrorCode());
			e.printStackTrace();
		}
	}

	private void initModel() {
		String err_msg = "Gurobi error when ";
		try {
			model = new GRBModel(env);
			err_msg = "creating the model. ";
			if (time_limit_min != null){
                model.set(GRB.DoubleParam.TimeLimit, time_limit_min * 60);
                err_msg += "setting the time limit. ";
			}

			if (node_limit != null){
				model.set(GRB.DoubleParam.NodeLimit, node_limit);
				err_msg += "setting the node limit. ";
            }

			if (max_solutions != null){
                model.set(GRB.IntParam.PoolSolutions, max_solutions);
                err_msg += "setting the maximum number of solutions. ";

				if (all_optimal) {
                    model.set(GRB.IntParam.PoolSearchMode, 2);
                    err_msg += "setting the pool search mode. ";
                }
            }
        } catch (GRBException e) {
			System.err.println(err_msg + " Error code: " + e.getErrorCode());
			e.printStackTrace();
		}
	}

	public void beginSearch() {
		super.beginSearch();
		initEnv();
		initModel();

		//create variable where solutions are stored
		GRBVar[] indicator_T_gur = new GRBVar[n];		
		for (int i = 0; i < n; i++) {
			try {
				indicator_T_gur[i] = model.addVar(0.0, 1.0, 0.0, GRB.BINARY, "w" + i);
			}
			catch (GRBException e) {
				System.err.println("Gurobi error when setting the time limit for observation #" + i + " Error code: " + e.getErrorCode());
				e.printStackTrace();
			}
		}
				
		//set objective
		GRBQuadExpr obj = new GRBQuadExpr();
		if (objective.equals(ObjectiveFunction.MAHAL)) {
			mahalDistSearch(model, indicator_T_gur, obj);
		}			
		else if (objective.equals(ObjectiveFunction.KER)) {
			kernelSearch(model, indicator_T_gur, obj);
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

		//ensures equal number of treatments and controls
		try {
			model.addConstr(expr, GRB.EQUAL, n / 2, "c0");
		} catch (GRBException e) {
			System.err.println("Gurobi error when setting the constraint of equal treatments and equal controls. Error code: " + e.getErrorCode());
			e.printStackTrace();
		}
		
		// Optimize model
		Thread optimize_thread = new Thread() {
			public void run() {
				try {
					model.optimize();
				} catch (GRBException e) {
					System.err.println("Gurobi error when running the optimization algorithm. Error code: " + e.getErrorCode());
					e.printStackTrace();
				}
			}
		};
		optimize_thread.start();
		try {
			optimize_thread.join();
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
         
        int num_solutions = 1;
        
        try {
			num_solutions = model.get(GRB.IntAttr.SolCount);
		} catch (GRBException e) {
			System.err.println("Gurobi error when querying the number of solutions. Error code: " + e.getErrorCode());
			e.printStackTrace();
		}
          
    	indicator_T = new int[num_solutions][n];
    	
        for (int k = 0; k < num_solutions; k++) {
        	try {
				model.set(GRB.IntParam.SolutionNumber, k);
			} catch (GRBException e) {
				System.err.println("Gurobi error when setting the solution number to: " + (k + 1) + " of " + num_solutions + " solutions. Error code: " + e.getErrorCode());
				e.printStackTrace();
			}
        	indicator_T_gur = model.getVars();
        	
            for (int i = 0; i < n; i++) {
        		indicator_T[k][i] = -99; //this is a "bad flag" to indicate to the user something went wrong
        	}
        	for (int i = 0; i < n; i++) {
        		try {
    				indicator_T[k][i] = (int)indicator_T_gur[i].get(GRB.DoubleAttr.X);
//        			System.out.println(indicator_T_gur[i].get(GRB.StringAttr.VarName) + " " +indicator_T_gur[i].get(GRB.DoubleAttr.X));
    			} catch (GRBException e) {
    				System.err.println("Gurobi error when extracting the solution for vector element #" + i + " for solution # " + (k + 1) + ". Error code: " + e.getErrorCode());
    				e.printStackTrace();
    			}
        	}
        }

		/*
		// for debugging- write model and solution to file
		try {
			model.write("model.lp");
			model.write("model.sol");
		} catch (GRBException e) {
			e.printStackTrace();
		}
		*/

		// Dispose of model and environment
        model.dispose();
        try {
			env.dispose();
		} catch (GRBException e) {
			System.err.println("Gurobi error when disposing of the environment. Error code: " + e.getErrorCode());
			e.printStackTrace();
		}
	}
	
	public void stopSearch() {
		super.stopSearch();
		model.terminate();
	}
	
	private void kernelSearch(GRBModel model, GRBVar[] indicator_T_gur, GRBQuadExpr obj) {
	
		//Setting objective matrix		
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				// 4 / n^2 * w^T K w = (2b - 1)^T K (2b - 1) = 
				// 4 / n^2 * (4 b^T K b - 4 1^T K b + 1^T K 1) \propto 
				// b^T K b - 1^T K b
				// sum_i sum_j K_ij b_i b_j - sum_i sum_j K_ij b_j
				obj.addTerm(Kgram[i][j], indicator_T_gur[i], indicator_T_gur[j]);	
				obj.addTerm(-Kgram[i][j], indicator_T_gur[j]);
			}
		}
	}

	private void mahalDistSearch(GRBModel model, GRBVar[] indicator_T_gur, GRBQuadExpr obj) {
    	
		SimpleMatrix Xsm = new SimpleMatrix(X);
		SimpleMatrix Sinvsm = new SimpleMatrix(Sinv);
		SimpleMatrix XSinvXt = Xsm.mult(Sinvsm).mult(Xsm.transpose());
		
		//Setting objective matrix
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				//hack to make number large here...
				obj.addTerm(1000000 * XSinvXt.get(i, j), indicator_T_gur[i], indicator_T_gur[j]);	    			
			}
		}
	}

	public void turnGurobiLogOff(){
		this.gurobi_log_off = true; 
	}
	
	public void setLogFilename(String log_file){
		this.log_file = log_file;
	}

	public void setTimeLimitMin(double time_limit_min) {
		this.time_limit_min = time_limit_min;
	}
	
	public void setTimeLimitMin(int node_limit) {
		this.node_limit = node_limit;
	}
	
	public void setMaxSolutions(int max_solutions){
		this.max_solutions = max_solutions;
	}
	
	public void setAllOptimal() {
		all_optimal = true;
	}
	
	public int[][] getIndicTs() {
		return indicator_T;
	}

}
