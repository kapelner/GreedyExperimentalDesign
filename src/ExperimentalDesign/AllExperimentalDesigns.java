package ExperimentalDesign;

import java.io.IOException;
import java.io.PrintStream;
import java.util.Random;
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

public abstract class AllExperimentalDesigns {

	//valid objective functions
	public static final String MAHAL = "mahal_dist";
	public static final String ABS = "abs_sum_diff";
	
	public Random r;
	
	//set by user
	protected int n;
	protected int p;
	protected String objective;
	protected Integer num_cores;
	
	//data inputed from the user's data
	protected double[][] Xstd;
	protected double[][] Sinv;	
	protected boolean wait;
	
	//temporary objects needed for search
	protected ExecutorService search_thread_pool;	
	protected boolean began_search;
	protected long t0;
	protected Long tf;
	
	public AllExperimentalDesigns(){
//		writeStdOutToLogFile();
//		System.out.println("GreedyExperimentalDesign");
		r = new Random();
	}	
	
	public void beginSearch(){
//		System.out.println("beginSearch");
		began_search = true;
		
		t0 = System.currentTimeMillis();
		//build the pool and all tasks to it
		search_thread_pool = Executors.newFixedThreadPool(num_cores == null ? 1 : num_cores);
	}
	
	
	protected void afterBeginSearch() {
		search_thread_pool.shutdown(); //run em all (but not on this thread!)
		Thread await_completion = new Thread(){
			public void run(){
				try {
					search_thread_pool.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS); //infinity
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
		search_thread_pool.shutdownNow();
	}	
	
	public void setObjective(String objective) throws Exception{
//		System.out.println("setObjective " + objective);
		this.objective = objective;
		if (!MAHAL.equals(objective) && !ABS.equals(objective)){
			throw new Exception("bad objective function");
		}
	}
	
	public void setNumCores(int num_cores){
//		System.out.println("setNumCores " +num_cores);
		this.num_cores = num_cores;
	}	
	
	public void setNandP(int n, int p) throws Exception{
		if (n % 2 != 0){
			throw new Exception("n must be even");
		}
//		System.out.println("setNandP n = " + n + " p = " + p);
		this.n = n;
		this.p = p;
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
	
	public void setSeed(int seed){
		r.setSeed(seed);
	}	
	
	public void setWait(){
		wait = true;
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
