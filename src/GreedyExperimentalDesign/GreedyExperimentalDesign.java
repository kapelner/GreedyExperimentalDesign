/*
    GreedyExperimentalDesign
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

package GreedyExperimentalDesign;

import java.util.ArrayList;
import ExperimentalDesign.AllExperimentalDesigns;
import ExperimentalDesign.Tools;

/**
 * This class handles initializing many greedy searches for a treatment vector
 * (the design) using a thread pool.
 * 
 * @author Adam Kapelner
 */
public class GreedyExperimentalDesign extends AllExperimentalDesigns {
	
	//set by user
	private int max_designs;
	private boolean diagnostics;
	private boolean semigreedy;

	//data inputed from the user's data
	private int[][] starting_indicTs;
	private Integer max_iters;
	
	//output
	private int[][] ending_indicTs;	
	private ArrayList<ArrayList<int[]>> switched_pairs;
	private ArrayList<ArrayList<double[]>> xbardiffjs_by_iterations;
	private Double[] objective_vals;
	private Integer[] num_iters;
	

	//running the Java as standalone is for debug purposes ONLY!!!
	public static void main(String[] args) throws Exception{	

		
		GreedyExperimentalDesign gd = new GreedyExperimentalDesign();
		//set seed here for reproducibility during debugging
		gd.r.setSeed(1984);

		int n = 1000;
		int p = 20;
		gd.setNandP(n, p);
		for (int i = 0; i < n; i++){
//			double[] x_i = {Math.random(), Math.random(), Math.random(), Math.random(), Math.random(), Math.random(), Math.random()};
			double[] x_i = new double[p];
			for (int j = 0; j < p; j++){
				x_i[j] = gd.r.nextDouble();
			}
			gd.setDataRow(i, x_i);
		}
//		System.out.println("Xstd");
//		for (int i = 0; i < n; i++){
//			System.out.println(Tools.StringJoin(gd.Xstd[i]));
//		}		
		gd.setMaxDesigns(25);
		gd.setObjective(ABS);
		gd.beginSearch();
//		System.out.println("progress: " + gd.progress());
	}
	
	public void beginSearch(){
		super.beginSearch();

		//initialize all data
		objective_vals = new Double[max_designs];
		num_iters = new Integer[max_designs];
		ending_indicTs = new int[max_designs][n];
		switched_pairs = new ArrayList<ArrayList<int[]>>(max_designs);
		for (int d = 0; d < max_designs; d++){
			switched_pairs.add(new ArrayList<int[]>());
		}
		xbardiffjs_by_iterations = new ArrayList<ArrayList<double[]>>(max_designs);
		for (int d = 0; d < max_designs; d++){
			xbardiffjs_by_iterations.add(new ArrayList<double[]>());
		}		
//		System.out.println("resulting data initialized");
		
		initializeStartingIndicTs();
		
		//convert Sinv to a matrix for easier multiplication inside the search
		final double[][] Sinvmat = new double[p][p];
		if (Sinv != null){
			for (int i = 0; i < p; i++){
				for (int j = 0; j < p; j++){
					Sinvmat[i][j] = Sinv[i][j];
				}			
			}
//			System.out.println("Sinvmat initialized");
		}
		

		for (int d = 0; d < max_designs; d++){
			final int d0 = d;
//			if (d % 100 == 0){
//				System.out.println("worker added to thread pool #" + d);
//			}
	    	search_thread_pool.execute(new Runnable(){
				public void run() {
					new GreedySearch(
							Xstd, 
							Sinvmat, 
							starting_indicTs[d0], 
							ending_indicTs[d0], 
							switched_pairs.get(d0),
							xbardiffjs_by_iterations.get(d0),
							objective_vals, 
							num_iters, 
							objective, 
							d0, 
							semigreedy, 
							diagnostics, 
							max_iters, 
							r);
				}
			});
		}		
		afterBeginSearch();
	}


	private void initializeStartingIndicTs() {
		starting_indicTs = new int[max_designs][n];
		for (int d = 0; d < max_designs; d++){
			starting_indicTs[d] = Tools.fisherYatesShuffle(Tools.newBlankDesign(n), r);
		}
	}

	public int progress(){
//		System.out.println("max_designs " + max_designs);
		int done = 0;
		if (objective_vals != null){
			for (int d = 0; d < max_designs; d++){
//				System.out.println("progress loop d = " + d);
				if (objective_vals[d] == null){
					break;
				}
				done++;
			}
		}
		return done;
	}
	
	public int[] getNumIters(){		
		int d_finished = progress();
		int[] num_iters = new int[d_finished];
		for (int i = 0; i < d_finished; i++){
			num_iters[i] = this.num_iters[i];
		}
		return num_iters;
	}
	
	public double[] getObjectiveVals(){		
		int d_finished = progress();
		double[] objective_vals = new double[d_finished];
		for (int i = 0; i < d_finished; i++){
			objective_vals[i] = this.objective_vals[i];
		}
		return objective_vals;
	}
	
	public int[][] getStartingIndicTs(int[] indicies){
		int[][] starting_indicTs = new int[indicies.length][n];
		for (int i = 0; i < indicies.length; i++){
			starting_indicTs[i] = this.starting_indicTs[indicies[i]];
		}
		return starting_indicTs;
	}
	
	public int[][] getEndingIndicTs(int[] indicies){
		int[][] ending_indicTs = new int[indicies.length][n];
		for (int i = 0; i < indicies.length; i++){
			ending_indicTs[i] = this.ending_indicTs[indicies[i]];
		}
		return ending_indicTs;
	}	
	
	public int[][][] getSwitchedPairs(int[] indicies){		
		int[][][] pairs_by_iteration_per_search = new int[indicies.length][][];
		for (int i = 0; i < indicies.length; i++){
			ArrayList<int[]> iteration_switched_pairs_raw = switched_pairs.get(indicies[i]);
			int iters = iteration_switched_pairs_raw.size();
			int[][] iteration_switched_pairs = new int[iters][];
			for (int j = 0; j < iters; j++){
				iteration_switched_pairs[j] = iteration_switched_pairs_raw.get(j);
			}
			pairs_by_iteration_per_search[i] = iteration_switched_pairs;
		}
		return pairs_by_iteration_per_search;
	}
	
	public double[][][] getXbarjDiffs(int[] indicies){		
		double[][][] xbarj_diffs_per_search = new double[indicies.length][][];
		for (int i = 0; i < indicies.length; i++){
			ArrayList<double[]> xbarj_diffs_raw = xbardiffjs_by_iterations.get(indicies[i]);
			int iters = xbarj_diffs_raw.size();
			double[][] xbarj_diffs = new double[iters][];
			for (int j = 0; j < iters; j++){
				xbarj_diffs[j] = xbarj_diffs_raw.get(j);
			}
			xbarj_diffs_per_search[i] = xbarj_diffs;
		}
		return xbarj_diffs_per_search;
	}	
		
	public void setMaxDesigns(int max_designs){
//		System.out.println("setMaxDesigns " + max_designs);
		this.max_designs = max_designs;
	}

	
	public void setDiagnostics(){
		diagnostics = true;
	}
	
	public void setSemigreedy(){
		semigreedy = true;
	}
	
	public void setMaxIters(int max_iters){
		this.max_iters = max_iters;
	}
}