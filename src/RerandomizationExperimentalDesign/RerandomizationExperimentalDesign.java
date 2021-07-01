package RerandomizationExperimentalDesign;

import java.util.ArrayList;

import ObjectiveFunctions.*;
import ExperimentalDesign.*;

public class RerandomizationExperimentalDesign extends MultipleSearchExperimentalDesigns {
	
	//set by user
	private Double obj_val_cutoff_to_include;
	
	//running the Java as standalone is for debug purposes ONLY!!!
	public static void main(String[] args) throws Exception{	

		RerandomizationExperimentalDesign rd = new RerandomizationExperimentalDesign();
		//set seed here for reproducibility during debugging
		rd.rand_obj.setSeed(1984);

		int n = 10;
		int p = 20;
		rd.setN(n);
		rd.setP(p);
		for (int i = 0; i < n; i++){
			double[] x_i = new double[p];
			for (int j = 0; j < p; j++){
				x_i[j] = rd.rand_obj.nextDouble();
			}
			rd.setDataRow(i, x_i);
		}
		
		rd.Sinv = new double[p][p];
		for (int i = 0; i < p; i++){
			for (int j = 0; j < p; j++){
				rd.Sinv[i][j] = i == j ? 1 : 0;			
			}
		}
		
		rd.setObjective(ObjectiveFunction.MAHAL);
		rd.setMaxDesigns(1000);
		rd.setObjValCutoffToInclude(10000D);
		rd.setNumCores(3);
		rd.setWait();
		rd.beginSearch();
		
		double[] obj_vals = rd.getObjectiveVals();
		System.out.println("obj_vals: " + Tools.StringJoin(obj_vals));	
		
		for (int i = 0; i < rd.ending_indicTs.length; i++){
			System.out.println("indicT " + (i + 1) + ": " + Tools.StringJoin(rd.ending_indicTs[i]));
		}
	
	}
		
	public void beginSearch(){
		super.beginSearch();
		
		//we gotta calculate the obj function
		ObjectiveFunction obj_fun = null;
		if (objective.equals(ObjectiveFunction.MAHAL)){
			obj_fun = new MahalObjective(Sinv, n);
		}
		else if (objective.equals(ObjectiveFunction.ABS)){
			obj_fun = new AbsSumObjective();	
		}
		else if (objective.equals(ObjectiveFunction.KER)){
			obj_fun = new KernelObjective(Kgram);	
		}
		
		final ObjectiveFunction fin_obj_fun = obj_fun;

//		System.out.println("before pool");

    	search_thread_pool.execute(new Runnable(){
			public void run() {
				while (true){
//					System.out.println("progress = " + r);
					//break up here too to avoid one more iteration (ugly, but a tad faster)
					if (num_completed.get() >= max_designs || search_stopped.get()){
						break;
					}
					
					int[] indicT = Tools.fisherYatesShuffle(Tools.newBalancedBlankDesign(n), rand_obj);
//					System.out.println("indicT " + Tools.StringJoin(indicT));
						
					if (objective.equals(ObjectiveFunction.KER)){
						((KernelObjective)fin_obj_fun).setW(indicT);							
					} else {
						int[] i_Ts = Tools.findIndicies(indicT, n / 2, 1);
						int[] i_Cs = Tools.findIndicies(indicT, n / 2, 0);
						ArrayList<double[]> XT = Tools.subsetMatrix(X, i_Ts); 
						ArrayList<double[]> XC = Tools.subsetMatrix(X, i_Cs);
						double[] avg_Ts = Tools.colAvg(XT, p);
						double[] avg_Cs = Tools.colAvg(XC, p);	
						fin_obj_fun.setXTbar(avg_Ts);
						fin_obj_fun.setXCbar(avg_Cs);
					}
					double obj_val = fin_obj_fun.calc(false);
					
					if (obj_val < obj_val_cutoff_to_include){
						synchronized(num_completed) {
							if (num_completed.get() >= max_designs || search_stopped.get()){
								break;
							}
							//create the new vector and its corresponding objective value
							ending_indicTs[num_completed.get()] = indicT;
							objective_vals[num_completed.get()] = obj_val;
							
//							System.out.println("num_completed: " + (num_completed + 1) + " obj val: " + obj_val);
//							System.out.println("w: " + Tools.StringJoin(indicT, ", "));
//							int[] i_Ts = Tools.findIndicies(indicT, n / 2, 1);
//							int[] i_Cs = Tools.findIndicies(indicT, n / 2, 0);
//							ArrayList<double[]> XT = Tools.subsetMatrix(X, i_Ts); 
//							ArrayList<double[]> XC = Tools.subsetMatrix(X, i_Cs);
//							
//							System.out.println("i_Ts size: " + i_Ts.length);
//							System.out.println("i_Cs size: " + i_Cs.length);
//							
//							System.out.println("XT size: " + XT.size());
//							System.out.println("XC size: " + XC.size());
//							
//							System.out.print("Ts: ");
//							double[] avg_Ts = Tools.colAvgDebug(XT, p);
//							System.out.print("\nCs: ");
//							double[] avg_Cs = Tools.colAvgDebug(XC, p);	
//							System.out.print("\n");
////							System.out.println("Ts: " + Tools.StringJoin(XT, ", ") + " Cs: " + Tools.StringJoin(XC, ", "));
//							System.out.println("xbarT: " + Tools.StringJoin(avg_Ts, " ") + " xbarC: " + Tools.StringJoin(avg_Cs, " "));
							num_completed.getAndIncrement();
						}
					}
				}
			}
    	});		
		afterBeginSearch();		
	}
	
	public void setObjValCutoffToInclude(double obj_val_cutoff_to_include){
		this.obj_val_cutoff_to_include = obj_val_cutoff_to_include;
	}	
}
