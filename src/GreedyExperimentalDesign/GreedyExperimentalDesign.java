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

import java.io.IOException;
import java.io.PrintStream;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.logging.FileHandler;
import java.util.logging.LogManager;
import java.util.logging.Logger;
import java.util.logging.StreamHandler;

import CustomLogging.LoggingOutputStream;
import CustomLogging.StdOutErrLevel;
import CustomLogging.SuperSimpleFormatter;
import no.uib.cipr.matrix.DenseMatrix;

/**
 * This class handles initializing many greedy searches for a treatment vector
 * (the design) using a thread pool.
 * 
 * @author Adam Kapelner
 */
public class GreedyExperimentalDesign {
	
	//valid objective functions
	public static final String MAHAL = "mahal_dist";
	public static final String ABS = "abs_sum_diff";
	
	//set by user
	private int n;
	private int p;
	private int max_designs;
	private boolean semigreedy;
	private String objective;
	private Integer num_cores;
	//data inputed from the user's data
	protected double[][] Xstd;
	private double[][] Sinv;
	private int[][] starting_indicTs;
	private boolean wait;
	private Integer max_iters;
	//temporary objects needed for search
	private ExecutorService greedy_search_thread_pool;
	private boolean began_search;
	private long t0;
	private Long tf;	
	//output
	private int[][] ending_indicTs;	
	private Double[] objective_vals;
	private Integer[] num_iters;

	//running the Java as standalone is for debug purposes ONLY!!!
	public static void main(String[] args) throws Exception{	
		
		GreedyExperimentalDesign gd = new GreedyExperimentalDesign();
		int n = 100;
		int p = 20;
		gd.setNandP(n, p);
		for (int i = 0; i < n; i++){
//			double[] x_i = {Math.random(), Math.random(), Math.random(), Math.random(), Math.random(), Math.random(), Math.random()};
			double[] x_i = new double[p];
			for (int j = 0; j < p; j++){
				x_i[j] = Math.random();
			}
			gd.setDataRow(i, x_i);
		}
//		System.out.println("Xstd");
//		for (int i = 0; i < n; i++){
//			System.out.println(Tools.StringJoin(gd.Xstd[i]));
//		}		
		gd.setMaxDesigns(100);
		gd.setObjective(MAHAL);
		gd.beginSearch();
		System.out.println("progress: " + gd.progress());
	}
	
	public GreedyExperimentalDesign(){
		writeStdOutToLogFile();
//		System.out.println("GreedyExperimentalDesign");
	}
	
	public void beginSearch(){
//		System.out.println("beginSearch");
		began_search = true;
		//initialize all data
		objective_vals = new Double[max_designs];
		num_iters = new Integer[max_designs];
		ending_indicTs = new int[max_designs][n];
//		System.out.println("resulting data initialized");
		
		initializeStartingIndicTs();
		
		//convert Sinv to a matrix for easier multiplication inside the search
		final DenseMatrix Sinvmat = new DenseMatrix(p, p);
		if (Sinv != null){
			for (int i = 0; i < p; i++){
				for (int j = 0; j < p; j++){
					Sinvmat.set(i, j, Sinv[i][j]);
				}			
			}
//			System.out.println("Sinvmat initialized");
		}
		
		t0 = System.currentTimeMillis();
		//build the pool and all tasks to it
		greedy_search_thread_pool = Executors.newFixedThreadPool(num_cores == null ? 1 : num_cores);
		for (int d = 0; d < max_designs; d++){
			final int d0 = d;
//			if (d % 100 == 0){
//				System.out.println("worker added to thread pool #" + d);
//			}
	    	greedy_search_thread_pool.execute(new Runnable(){
				public void run() {
					new GreedySearch(Xstd, Sinvmat, starting_indicTs[d0], ending_indicTs[d0], objective_vals, num_iters, objective, d0, semigreedy, max_iters);
				}
			});
		}
		greedy_search_thread_pool.shutdown(); //run em all (but not on this thread!)
		Thread await_completion = new Thread(){
			public void run(){
				try {
					greedy_search_thread_pool.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS); //infinity
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
				tf = System.currentTimeMillis();
			}
		};
		await_completion.start();
		if (wait){
			try {
				await_completion.join();
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
	}
	
	private void initializeStartingIndicTs() {
		starting_indicTs = new int[max_designs][n];
		for (int d = 0; d < max_designs; d++){
			starting_indicTs[d] = Tools.fisherYatesShuffle(Tools.newBlankDesign(n));
		}
	}

	public int timeElapsedInSeconds(){
		if (tf == null){
			return (int)(System.currentTimeMillis() - t0) / 1000;
		}
		return (int)(tf - t0) / 1000;
	}
	
	public long timeFinished(){
		return tf;
	}
	
	public boolean began(){
		return began_search;
	}
	
	public void stopSearch(){
		greedy_search_thread_pool.shutdownNow();
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
	
	public int[][] getEndingIndicTs(int[] indicies){
		int[][] ending_indicTs = new int[indicies.length][n];
		for (int i = 0; i < indicies.length; i++){
			ending_indicTs[i] = this.ending_indicTs[indicies[i]];
		}
		return ending_indicTs;
	}	
	
	public void setMaxDesigns(int max_designs){
//		System.out.println("setMaxDesigns " + max_designs);
		this.max_designs = max_designs;
	}
	
	public void setNumCores(int num_cores){
//		System.out.println("setNumCores " +num_cores);
		this.num_cores = num_cores;
	}	
	
	public void setNandP(int n, int p){
//		System.out.println("setNandP n = " + n + " p = " + p);
		this.n = n;
		this.p = p;
	}
	
	public void setObjective(String objective) throws Exception{
//		System.out.println("setObjective " + objective);
		this.objective = objective;
		if (!MAHAL.equals(objective) && !ABS.equals(objective)){
			throw new Exception("bad objective function");
		}
	}
	
	public void setDataRow(int i0, double[] x_i){
//		System.out.println("setDataRow " + i0 + "  " + x_i);
		if (Xstd == null){
			Xstd = new double[n][p];
		}
		for (int j = 0; j < p; j++){
			Xstd[i0][j] = x_i[j];
		}
	}
	
	public void setDataRow(int i0, double x_i){
//		System.out.println("setDataRow " + i0 + "  " + x_i);
		if (Xstd == null){
			Xstd = new double[n][p];
		}
		double[] row = {x_i};
		Xstd[i0] = row;
	}	
	
	public void setInvVarCovRow(int j0, double[] Sinv_i){
//		System.out.println("setInvVarCovRow " + j0 + "  " + Sinv_i);
		if (Sinv == null){
			Sinv = new double[p][p];
		}
		for (int j = 0; j < p; j++){
			Sinv[j0][j] = Sinv_i[j];
		}
	}
	
	public void setInvVarCovRow(int j0, double Sinv_i){
//		System.out.println("setInvVarCovRow " + j0 + "  " + Sinv_i);
		if (Sinv == null){
			Sinv = new double[p][p];
		}
		double[] row = {Sinv_i};
		Sinv[j0] = row;
	}
	
	public void setSemigreedy(){
		semigreedy = true;
	}
	
	public void setWait(){
		wait = true;
	}
	
	public void setMaxIters(int max_iters){
		this.max_iters = max_iters;
	}
	
	public static void writeStdOutToLogFile(){
		try {
		  Logger.getLogger("").addHandler(new StreamHandler()); //turn off std out
		  suppressOrWriteToDebugLog();
		}
		catch (Error e){
			System.out.println("Logger and or suppressOrWriteToDebugLog FAILING\n");
		}    
 	}
	
	public static void suppressOrWriteToDebugLog(){
		//also handle the logging
        LogManager logManager = LogManager.getLogManager();
        logManager.reset();

        // create log file, no limit on size
        FileHandler fileHandler = null;
		try {
			fileHandler = new FileHandler("java_log" + ".log", Integer.MAX_VALUE, 1, false);
		} catch (SecurityException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
        fileHandler.setFormatter(new SuperSimpleFormatter());
        Logger.getLogger("").addHandler(fileHandler);
        
        
        // now rebind stdout/stderr to logger
        Logger logger = Logger.getLogger("stdout");         
        LoggingOutputStream  los = new LoggingOutputStream(logger, StdOutErrLevel.STDOUT);
        System.setOut(new PrintStream(los, true));
        logger = Logger.getLogger("stderr");                                    
        los = new LoggingOutputStream(logger, StdOutErrLevel.STDERR);            
        System.setErr(new PrintStream(los, true)); 		
	}
}