/*
    OptimalExperimentalDesign
    Software for Experimental Design
    
    Copyright (C) 2015 Professor Adam Kapelner 
    Department of Mathematics, Queens College, City University of New York

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details:
    
    http://www.gnu.org/licenses/gpl-2.0.txt

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

package OptimalExperimentalDesign;

import ExperimentalDesign.AbsSumObjective;
import java.util.ArrayList;
import ExperimentalDesign.AllExperimentalDesigns;
import ExperimentalDesign.ObjectiveFunction;
import ExperimentalDesign.PropMahalObjective;
import ExperimentalDesign.Tools;
import GreedyExperimentalDesign.GreedyExperimentalDesign;

/**
 * This class handles initializing many optimal searches for a treatment vector
 * (the design) using a thread pool.
 * 
 * @author Adam Kapelner
 */
public class OptimalExperimentalDesign extends AllExperimentalDesigns {

	private static final double NOT_REACHED_YET = -999999;
	
	//valid objective functions
	public static final String MAHAL = "mahal_dist";
	public static final String ABS = "abs_sum_diff";
	
	//temp stuff
	private ArrayList<int[]> all_indicTs;
	private int max_designs;
	private int n_over_two;
	
	//output
	private double[] objective_vals;
	private int d_opt;

	//running the Java as standalone is for debug purposes ONLY!!!
	public static void main(String[] args) throws Exception{	

		OptimalExperimentalDesign od = new OptimalExperimentalDesign();
		//set seed here for reproducibility during debugging
		od.r.setSeed(1984);

		int n = 28;
		int p = 1;
		od.setNandP(n, p);
		for (int i = 0; i < n; i++){
//			double[] x_i = {Math.random(), Math.random(), Math.random(), Math.random(), Math.random(), Math.random(), Math.random()};
			double[] x_i = new double[p];
			for (int j = 0; j < p; j++){
				x_i[j] = od.r.nextDouble();
			}
			od.setDataRow(i, x_i);
		}
//		System.out.println("Xstd");
//		for (int i = 0; i < n; i++){
//			System.out.println(Tools.StringJoin(od.Xstd[i]));
//		}
		od.setObjective(ABS);
		od.setNumCores(2);
		od.beginSearch();
		od.setWait();
//		System.out.println("progress: " + od.progress());
	}
	
	public void beginSearch(){
		super.beginSearch();
		
		initializeStartingIndicTs();
		
		objective_vals = new double[max_designs];
		for (int d = 0; d < max_designs; d++){
			objective_vals[d] = NOT_REACHED_YET;
		}
		
		for (int d = 0; d < max_designs; d++){
			final int d0 = d;
//			if (d % 100 == 0){
//				System.out.println("worker added to thread pool #" + d);
//			}
	    	search_thread_pool.execute(new Runnable(){
				public void run() {
					ObjectiveFunction obj_fun = null;
					if (objective.equals(GreedyExperimentalDesign.MAHAL)){
						obj_fun = new PropMahalObjective(Sinv);
					}
					else if (objective.equals(GreedyExperimentalDesign.ABS)){
						obj_fun = new AbsSumObjective();	
					}
					
					//get the vector for this run
					int[] indicT = all_indicTs.get(d0);
					//get the indicies
					int[] i_Ts = Tools.findIndicies(indicT, n_over_two, 1);
//					System.out.println("i_Ts " + Tools.StringJoin(i_Ts));
					int[] i_Cs = Tools.findIndicies(indicT, n - n_over_two, 0);
					//get the rows for each group
					ArrayList<double[]> XT = Tools.subsetMatrix(Xstd, i_Ts); 
					ArrayList<double[]> XC = Tools.subsetMatrix(Xstd, i_Cs); 
					//compute the averages
					double[] avg_Ts = Tools.colAvg(XT, p);
					double[] avg_Cs = Tools.colAvg(XC, p);
					obj_fun.setXTbar(avg_Ts);
					obj_fun.setXCbar(avg_Cs);
					
					//calculate our objective function (according to the user's specification)
					objective_vals[d0] = obj_fun.calc(false);
				}
			});
		}
		afterBeginSearch();
		
//		for (int i = 0; i < max_designs; i++){
//			System.out.println((i + 1) + ": " + objective_vals[i]);
//		}
		//now return the min
		d_opt = Tools.min_index(objective_vals);
		System.out.println("d_opt: " + (d_opt + 1));
		System.out.println("obj_opt: " + objective_vals[d_opt]);
		System.out.println("time elapsed in sec: " + timeElapsedInSeconds());
	}
	
	private void initializeStartingIndicTs() {
		System.out.println("begin initializeStartingIndicTs");
		max_designs = (int)n_choose_k(n, n / 2);
		all_indicTs = new ArrayList<int[]>(max_designs);
//		int[] start = new int[n];
//		for (int i = 0; i < n; i++){
//			start[i] = 3;
//		}
//		recursivelyFindAllBinaryVecs(start, 0, 0, 0);
//		for (int i = 0; i < max_designs; i++){
//			System.out.println((i + 1) + ": " + Tools.StringJoin(all_indicTs.get(i), ""));
//		}
		recursivelyFindAllBinaryVecs(new int[n], 0, 0, 0);
		System.out.println("end initializeStartingIndicTs");
	}

	private void recursivelyFindAllBinaryVecs(int[] vec, int pos, int on, int off) {
		//if we've made it to the end, we're done
		if (pos == n){
			all_indicTs.add(vec);
			return;
		}
		
		//now set the next position on and recurse
		int[] next_vec_on = vec.clone();
		next_vec_on[pos] = 1;
		if (on + 1 <= n_over_two){
			recursivelyFindAllBinaryVecs(next_vec_on, pos + 1, on + 1, off);
		}		
		
		//now set the next position off and recurse
		int[] next_vec_off = vec.clone();
		next_vec_off[pos] = 0;
		if (off + 1 <= n_over_two){
			recursivelyFindAllBinaryVecs(next_vec_off, pos + 1, on, off + 1);
		}
	}

	private long n_choose_k(int n, int k) {		
		 long n_choose_k = (long)Math.exp(ln_factorial(n) - ln_factorial(k) - ln_factorial(n - k));
//		 System.out.println(n + " choose " + k + " = " + n_choose_k);
		 return n_choose_k;
	}

	private double ln_factorial(int n) {
		double sum = 0;
		for (int i = 1; i <= n; i++){
			sum += Math.log(i);
		}
//		System.out.println("ln " + n + "! = " + sum);
		return sum;
	}

	public int progress(){
//		System.out.println("max_designs " + max_designs);
		int done = 0;
		if (objective_vals != null){
			for (int d = 0; d < max_designs; d++){
//				System.out.println("progress loop d = " + d);
				if (objective_vals[d] == NOT_REACHED_YET){
					break;
				}
				done++;
			}
		}
		return done;
	}
	
	public double[] getObjectiveVals(){		
		int d_finished = progress();
		double[] objective_vals = new double[d_finished];
		for (int i = 0; i < d_finished; i++){
			objective_vals[i] = this.objective_vals[i];
		}
		return objective_vals;
	}
	
	
	public void setNandP(int n, int p) throws Exception {
		super.setNandP(n, p);
		n_over_two = n / 2;
	}

	

}