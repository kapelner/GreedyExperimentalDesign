package GurobiNumericalOptimizeExperimentalDesign;

import ExperimentalDesign.AllExperimentalDesigns;
import gurobi.*;

public class GurobiNumericalOptimizeExperimentalDesign extends AllExperimentalDesigns {

	//running the Java as standalone is for debug purposes ONLY!!!
	public static void main(String[] args) throws Exception{	

		GurobiNumericalOptimizeExperimentalDesign gnuoed = new GurobiNumericalOptimizeExperimentalDesign();
		//set seed here for reproducibility during debugging
		gnuoed.rand_obj.setSeed(1984);

		int n = 6;
		int p = 3;
		gnuoed.setNandP(n, p);
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
	
	public void beginSearch(){
		super.beginSearch();
		
		//bracha: the data "X" is the data matrix and "Sinv" is the inverse sample var-cov matrix"
				
	    try {
	        GRBEnv    env   = new GRBEnv("gurobi_numerical_optimization_via_R_package_Gree.log");
	        GRBModel  model = new GRBModel(env);

	        // Create variables

	        GRBVar x = model.addVar(0.0, 1.0, 0.0, GRB.BINARY, "x");
	        GRBVar y = model.addVar(0.0, 1.0, 0.0, GRB.BINARY, "y");
	        GRBVar z = model.addVar(0.0, 1.0, 0.0, GRB.BINARY, "z");

	        // Set objective: maximize x + y + 2 z

	        GRBLinExpr expr = new GRBLinExpr();
	        expr.addTerm(1.0, x); expr.addTerm(1.0, y); expr.addTerm(2.0, z);
	        model.setObjective(expr, GRB.MAXIMIZE);

	        // Add constraint: x + 2 y + 3 z <= 4

	        expr = new GRBLinExpr();
	        expr.addTerm(1.0, x); expr.addTerm(2.0, y); expr.addTerm(3.0, z);
	        model.addConstr(expr, GRB.LESS_EQUAL, 4.0, "c0");

	        // Add constraint: x + y >= 1

	        expr = new GRBLinExpr();
	        expr.addTerm(1.0, x); expr.addTerm(1.0, y);
	        model.addConstr(expr, GRB.GREATER_EQUAL, 1.0, "c1");

	        // Optimize model

	        model.optimize();

	        System.out.println(x.get(GRB.StringAttr.VarName)
	                           + " " +x.get(GRB.DoubleAttr.X));
	        System.out.println(y.get(GRB.StringAttr.VarName)
	                           + " " +y.get(GRB.DoubleAttr.X));
	        System.out.println(z.get(GRB.StringAttr.VarName)
	                           + " " +z.get(GRB.DoubleAttr.X));

	        System.out.println("Obj: " + model.get(GRB.DoubleAttr.ObjVal));

	        // Dispose of model and environment

	        model.dispose();
	        env.dispose();

	      } catch (GRBException e) {
	        System.out.println("Error code: " + e.getErrorCode() + ". " +
	                           e.getMessage());
	      }
	}
}
