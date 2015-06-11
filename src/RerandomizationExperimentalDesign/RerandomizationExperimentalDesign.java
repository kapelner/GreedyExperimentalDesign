package RerandomizationExperimentalDesign;

import java.util.ArrayList;

import ExperimentalDesign.AbsSumObjective;
import ExperimentalDesign.AllExperimentalDesigns;
import ExperimentalDesign.MahalObjective;
import ExperimentalDesign.ObjectiveFunction;
import ExperimentalDesign.Tools;

public class RerandomizationExperimentalDesign extends AllExperimentalDesigns {
	
	//set by user
	private int max_designs;
	private Double obj_val_cutoff_to_include;
	private ArrayList<IndicTAndObjVal> ending_indicTs_and_obj_vals;
	
	private class IndicTAndObjVal {
		public int[] indicT;
		public Double obj_val;
		
		public IndicTAndObjVal(int[] indicT, Double obj_val){
			this.indicT = indicT;
			this.obj_val = obj_val;
		}
	}
	
	//running the Java as standalone is for debug purposes ONLY!!!
	public static void main(String[] args) throws Exception{	

		RerandomizationExperimentalDesign rd = new RerandomizationExperimentalDesign();
		//set seed here for reproducibility during debugging
		rd.rand_obj.setSeed(1984);

		int n = 100;
		int p = 20;
		rd.setNandP(n, p);
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
		
		rd.setObjective(MAHAL);
		rd.setMaxDesigns(1000);
		rd.setObjValCutoffToInclude(0.3);
		rd.setNumCores(3);
		rd.setWait();
		rd.beginSearch();
		
		double[] obj_vals = rd.getObjectiveVals();
		System.out.println("obj_vals: " + Tools.StringJoin(obj_vals));	
		
		int[][] ending_indicTs = rd.getEndingIndicTs();	
		for (int i = 0; i < ending_indicTs.length; i++){
			System.out.println("indicT " + (i + 1) + ": " + Tools.StringJoin(ending_indicTs[i]));
		}
	
	}
		
	public void beginSearch(){
		super.beginSearch();
		
		ending_indicTs_and_obj_vals = new ArrayList<IndicTAndObjVal>(max_designs);
		

    	search_thread_pool.execute(new Runnable(){
			public void run() {
				while (true){
					//break up here too to avoid one more iteration (ugly, but a tad faster)
					if (ending_indicTs_and_obj_vals.size() == max_designs){
						break;
					}
					
					int[] indicT = Tools.fisherYatesShuffle(Tools.newBalancedBlankDesign(n), rand_obj);
					
					if (obj_val_cutoff_to_include != null){
						//we gotta calculate the obj function
						ObjectiveFunction obj_fun = null;
						if (objective.equals(MAHAL)){
							obj_fun = new MahalObjective(Sinv, n);
						}
						else if (objective.equals(ABS)){
							obj_fun = new AbsSumObjective();	
						}
						int[] i_Ts = Tools.findIndicies(indicT, n / 2, 1);
						int[] i_Cs = Tools.findIndicies(indicT, n / 2, 0);
						ArrayList<double[]> XT = Tools.subsetMatrix(X, i_Ts); 
						ArrayList<double[]> XC = Tools.subsetMatrix(X, i_Cs);
						double[] avg_Ts = Tools.colAvg(XT, p);
						double[] avg_Cs = Tools.colAvg(XC, p);	
						obj_fun.setXTbar(avg_Ts);
						obj_fun.setXCbar(avg_Cs);
						double obj_val = obj_fun.calc(false);
						
						if (obj_val < obj_val_cutoff_to_include){
							if (ending_indicTs_and_obj_vals.size() == max_designs){
								break;
							}
							ending_indicTs_and_obj_vals.add(new IndicTAndObjVal(indicT, obj_val));
//								if (ending_indicTs_and_obj_vals.size() % 10 == 0){
//									System.out.println(ending_indicTs_and_obj_vals.size() + " vectors found!");
//								}
						}
					}
					else {
						//we are just looking for a certain number and then we're done
						if (ending_indicTs_and_obj_vals.size() == max_designs){
							break;
						}
						ending_indicTs_and_obj_vals.add(new IndicTAndObjVal(indicT, null)); //we don't care about this
					}
					//break out if user desires
					if (search_thread_pool.isShutdown()){
						break;
					}
				}
			}
		});		
		afterBeginSearch();		
	}
	
	public void setMaxDesigns(int max_designs){
		this.max_designs = max_designs;
	}
	
	public void setObjValCutoffToInclude(double obj_val_cutoff_to_include){
		this.obj_val_cutoff_to_include = obj_val_cutoff_to_include;
	}
	
	public int progress(){
		return ending_indicTs_and_obj_vals.size();
	}
	
	
	public double[] getObjectiveVals(){
		int r = progress();
		double[] objective_vals = new double[r];
		for (int i = 0; i < r; i++){
			objective_vals[i] = ending_indicTs_and_obj_vals.get(i).obj_val;
		}
		return objective_vals;
	}
	
	public int[][] getEndingIndicTs(){
		int r = progress();
		int[][] ending_indicTs = new int[r][n];
		for (int i = 0; i < r; i++){
			ending_indicTs[i] = ending_indicTs_and_obj_vals.get(i).indicT;
		}
		return ending_indicTs;
	}	
}
