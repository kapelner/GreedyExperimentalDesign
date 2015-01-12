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
 * This class handles the parallelization of many Gibbs chains over many CPU cores
 * to create one BART regression model. It also handles all operations on the completed model.
 * 
 * @author Adam Kapelner and Justin Bleich
 */
public class GreedyExperimentalDesign {
	
	//valid objective functions
	public static final String MAHAL = "mahal_dist";
	public static final String ABS = "abs_sum_diff";
	
	//set by user
	private int n;
	private int p;
	private int max_designs;
	private String objective;
	private int num_cores;
	//data inputted from the user's data
	private double[][] Xstd;
	private double[][] Sinv;
	private int[][] starting_indicTs;
	
	//temp objects needed for search
	private ExecutorService greedy_search_thread_pool;
	
	//output
	private int[][] ending_indicTs;	
	private Double[] objective_vals;
	
	public GreedyExperimentalDesign(){
		System.out.println("GreedyExperimentalDesign");
//		writeStdOutToLogFile();
	}
	
	public void beginSearch(){
		//initialize all data
		objective_vals = new Double[max_designs];
		ending_indicTs = new int[max_designs][n];
		
		//convert Sinv to a matrix for easier multiplication inside the search
		final DenseMatrix Sinvmat = new DenseMatrix(p, p);
		for (int i = 0; i < p; i++){
			for (int j = 0; j < p; j++){
				Sinvmat.set(i, j, Sinv[i][j]);
			}			
		}
		
		//build the pool and all tasks to it
		greedy_search_thread_pool = Executors.newFixedThreadPool(num_cores);
		for (int d = 0; d < max_designs; d++){
			final int d0 = d;
	    	greedy_search_thread_pool.execute(new Runnable(){
				public void run() {
					new GreedySearch(Xstd, Sinvmat, starting_indicTs[d0], ending_indicTs[d0], objective_vals, objective, d0);
				}
			});
		}
		greedy_search_thread_pool.shutdown();
		try {	         
	         greedy_search_thread_pool.awaitTermination(Long.MAX_VALUE, TimeUnit.SECONDS); //effectively infinity
	    } catch (InterruptedException ignored){}	
		
	}
	
	public void stopSearch(){
		greedy_search_thread_pool.shutdownNow();
	}
	

	public int progress(){
		int done = 0;
		if (objective_vals != null){
			for (int d = 0; d < max_designs; d++){
				if (objective_vals[d] != null){
					done++;
				}
			}
		}
		return done;
	}
	
	public Double[] getObjectiveVals(){
		return objective_vals;
	}
	
	public int[] getEndingIndicT(int d){
		return ending_indicTs[d];
	}
	
	public void setMaxDesigns(int max_designs){
		System.out.println("setMaxDesigns " + max_designs);
		this.max_designs = max_designs;
	}
	
	public void setNumCores(int num_cores){
		System.out.println("setNumCores " +num_cores);
		this.num_cores = num_cores;
	}	
	
	public void setNandP(int n, int p){
		System.out.println("setNandP n = " + n + " p = " + p);
		this.n = n;
		this.p = p;
	}
	
	public void setObjective(String objective) throws Exception{
		System.out.println("setObjective " + objective);
		this.objective = objective;
		if (!MAHAL.equals(objective) && !ABS.equals(objective)){
			throw new Exception("bad objective function");
		}
	}
	
	public void setDataRow(int i0, double[] x_i){
		System.out.println("setDataRow " + i0 + "  " + x_i);
		if (Xstd == null){
			Xstd = new double[n][p];
		}
		for (int j = 0; j < p; j++){
			Xstd[i0][j] = x_i[j];
		}
	}
	
	public void setInvVarCovRow(int j0, double[] Sinv_i){
		System.out.println("setInvVarCovRow " + j0 + "  " + Sinv_i);
		if (Sinv == null){
			Sinv = new double[p][p];
		}
		for (int j = 0; j < p; j++){
			Sinv[j0][j] = Sinv_i[j];
		}
	}
	
	public void setDesignStartingPoint(int d0, int[] indicT){
		System.out.println("setDesignStartingPoint " + d0 + " " + indicT);
		if (starting_indicTs == null){
			starting_indicTs = new int[max_designs][n];
		}
		for (int i = 0; i < n; i++){
			starting_indicTs[d0][i] = indicT[i];
		}		
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