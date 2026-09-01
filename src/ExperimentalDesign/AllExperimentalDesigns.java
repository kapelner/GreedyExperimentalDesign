package ExperimentalDesign;

import java.util.Random;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicReference;

import CustomLogging.FileLoggedClass;
import ObjectiveFunctions.ObjectiveFunction;

public abstract class AllExperimentalDesigns extends FileLoggedClass {
	
	public Random rand_obj;
	
	//set by user
	protected int n;
	protected int p;
	protected String objective;
	protected Integer num_cores;
	protected AtomicBoolean search_stopped;
	protected AtomicReference<Throwable> worker_error;
	
	//data inputed from the user's data
	protected double[][] X;
	protected double[][] Sinv;	
	protected double[][] Kgram;
	protected boolean wait;
	protected boolean verbose = true;
	protected boolean use_gpu = false;
	
	//temporary objects needed for search
	protected ExecutorService search_thread_pool;	
	protected boolean began_search;
	protected long t0;
	protected Long tf;
	
	
	public AllExperimentalDesigns(){
		num_cores = 1;
		rand_obj = new Random();
		search_stopped = new AtomicBoolean();
		worker_error = new AtomicReference<>();
	}	
	
	public void beginSearch() {
//		System.out.println("beginSearch");
		began_search = true;
		worker_error.set(null);
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
			Throwable err = worker_error.get();
			if (err != null) {
				throw new RuntimeException("Search worker failed: " + err.getMessage(), err);
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
		search_stopped.set(true);
	}	
	
	public void setObjective(String objective) throws Exception{
		if (!ObjectiveFunction.isValidObjFunction(objective)){
			throw new Exception("objective function not recognized");
		}
		this.objective = objective;
	}
	
	public void setNumCores(int num_cores){
//		System.out.println("setNumCores " +num_cores);
		this.num_cores = num_cores;
	}	

	public void setVerbose(boolean verbose){
		this.verbose = verbose;
	}

	public boolean isVerbose(){
		return verbose;
	}

	public void setUseGpu(boolean use_gpu){
		this.use_gpu = use_gpu;
	}

	public boolean isUseGpu(){
		return use_gpu;
	}
	
	public void setN(int n) {
		this.n = n;
	}	
		
	public void setP(int p){
		this.p = p;
	}	
	
	public void setKgramRow(int i0, double[] kgram_i){
		if (Kgram == null){
			Kgram = new double[n][n];
		}
		for (int j = 0; j < n; j++){
			Kgram[i0][j] = kgram_i[j];
		}
//		System.out.println("setKgramRow " + i0 + " " + Tools.StringJoin(Kgram[i0]));
	}

	public void setKgramMatrix(double[] kgram_flat, int n){
		this.n = n;
		Kgram = new double[n][n];
		int idx = 0;
		for (int j = 0; j < n; j++){
			for (int i = 0; i < n; i++){
				Kgram[i][j] = kgram_flat[idx++];
			}
		}
	}
	
	public void setDataRow(int i0, double[] x_i){
//		System.out.println("setDataRow " + i0 + "  " + x_i);
		if (X == null){
			X = new double[n][p];
		}
		for (int j = 0; j < p; j++){
			X[i0][j] = x_i[j];
		}
	}

	public void setDataMatrix(double[] x_flat, int n, int p){
		this.n = n;
		this.p = p;
		X = new double[n][p];
		int idx = 0;
		for (int j = 0; j < p; j++){
			for (int i = 0; i < n; i++){
				X[i][j] = x_flat[idx++];
			}
		}
	}
	
	public void setDataRow(int i0, double x_i){
//		System.out.println("setDataRow " + i0 + "  " + x_i);
		if (X == null){
			X = new double[n][p];
		}
		double[] row = {x_i};
		X[i0] = row;
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

	public void setInvVarCovMatrix(double[] sinv_flat, int p){
		this.p = p;
		Sinv = new double[p][p];
		int idx = 0;
		for (int j = 0; j < p; j++){
			for (int i = 0; i < p; i++){
				Sinv[i][j] = sinv_flat[idx++];
			}
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
		rand_obj.setSeed(seed);
	}	
	
	public void setWait(){
		wait = true;
	}	
}
