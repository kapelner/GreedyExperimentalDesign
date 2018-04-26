package GurobiNumericalOptimizeExperimentalDesign;

import ExperimentalDesign.AllExperimentalDesigns;
import org.ejml.simple.SimpleMatrix;
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
	        
	        int n = X.length;
	    	int p = X[1].length;
	    	
	    	SimpleMatrix Xsm = new SimpleMatrix(X);
	    	SimpleMatrix Sinvsm = new SimpleMatrix(p, p);
	    	// Means:
	        SimpleMatrix x = new SimpleMatrix(p, 1);
	        for(int r=0; r<p; r++ ){
	            x.set(r, 0, Xsm.transpose().extractVector(true, r).elementSum() / n);
	        }
	        // Covariance matrix:
	        for(int r=0; r<p; r++){
	            for(int c=0; c<p; c++){
	                if(r > c){
	                    Sinvsm.set(r, c, Sinvsm.get(c, r));
	                } else {
	                    double cov = Xsm.transpose().extractVector(true, r).minus( x.get((r), 0) ).dot(Xsm.transpose().extractVector(true, c).minus( x.get((c), 0) ).transpose());
	                    Sinvsm.set(r, c, (cov / n));
	                }
	            }
	        }
	    	SimpleMatrix XSinvXt = Xsm.mult(Sinvsm).mult(Xsm.transpose());
	    	
	    	System.out.println("sinv: " + Sinvsm.get(0,0) + " n: "+n+" p : "+p);
	    	

	        // Create variables
	    	GRBVar[] vars = new GRBVar[n];
	    	for(int i = 0; i<n; i++) {
	    		vars[i] = model.addVar(0.0, 2.0, 0.0, GRB.BINARY, Integer.toString(i));
	    	}
	    	
	    	//Setting objective matrix
	    	GRBQuadExpr obj = new GRBQuadExpr();
	    	for(int i = 0; i<n; i++) {
	    		for(int j = 0; j<n; j++) {
	    			obj.addTerm(XSinvXt.get(i ,j), vars[i], vars[j]);	    			
	    		}
	    	}
	    	model.setObjective(obj);


	        // Add constraint: sum of vars equal to n/2 (equal num 1's and 0's)
	    	GRBLinExpr expr = new GRBLinExpr();
	    	for(int i = 0; i < n; i++) {
	    		 expr.addTerm(1.0, vars[i]);
	    	}
	        model.addConstr(expr, GRB.EQUAL, (n/2), "c0");

	        // Optimize model
	        model.optimize();
	        
	        for(int i = 0; i < n; i++) {
	        	System.out.println(vars[i].get(GRB.StringAttr.VarName)
                        + " " +vars[i].get(GRB.DoubleAttr.X));
	    	}

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
