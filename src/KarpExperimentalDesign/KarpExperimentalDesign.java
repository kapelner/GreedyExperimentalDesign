package KarpExperimentalDesign;

import ExperimentalDesign.AllExperimentalDesigns;
import ExperimentalDesign.Tools;

public class KarpExperimentalDesign extends AllExperimentalDesigns {
	
	
	private boolean balanced;
	private KarpDesignSearcher keds;
	
	//running the Java as standalone is for debug purposes ONLY!!!
	public static void main(String[] args) throws Exception{	

		for (int g = 0; g < 1; g++){
			KarpExperimentalDesign kd = new KarpExperimentalDesign();
			//set seed here for reproducibility during debugging
			kd.r.setSeed(1984);
	
			int n = 82;
			kd.setNandP(n, 1);
			for (int i = 0; i < n; i++){
				double[] x_i = new double[1];
				x_i[0] = kd.r.nextDouble();
//				x_i[0] = n - (i);
				kd.setDataRow(i, x_i);
			}
	//		System.out.println("Xstd");
	//		for (int i = 0; i < n; i++){
	//			System.out.println(Tools.StringJoin(od.Xstd[i]));
	//		}
			kd.setObjective(ABS);
			kd.setWait();
			kd.setBalanced();
			kd.beginSearch();
		}
//		System.out.println("progress: " + od.progress());
	}	
	
	public void beginSearch(){
		super.beginSearch();
    	search_thread_pool.execute(new Runnable(){
			public void run() {
				keds = balanced ? new KarpDesignSearcherBalanced(Xstd) : new KarpDesignSearcherUnbalanced(Xstd);
			}
		});
		afterBeginSearch();	
//		System.out.println("FINAL INDIC_T: " + Tools.StringJoin(getKarpIndicT()));
//		System.out.println("Num T: " + Tools.sum_array(getKarpIndicT()) + " n: " + n);
//		System.out.println("Final obj val: " + getKarpObjectiveVal());
	}
	
	public void setBalanced(){
		balanced = true;
	}
	
	public double getKarpObjectiveVal(){		
		return keds.getObjVal();
	}
	
	public int[] getKarpIndicT() {
		return keds.getIndicT();
	}	
	
	public double progress(){
		return keds == null ? 0 : keds.progress();
	}	
}
