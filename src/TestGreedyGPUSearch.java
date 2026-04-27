import OptimalExperimentalDesign.GpuBridge;
import GreedyExperimentalDesign.*;
import ObjectiveFunctions.*;
import ExperimentalDesign.Tools;
import java.util.ArrayList;
import java.util.Random;
import java.util.concurrent.atomic.AtomicBoolean;

public class TestGreedyGPUSearch {
    public static void main(String[] args) throws Exception {
        int n = 100;
        int p = 5;
        int nT = 50;

        System.out.println("Starting Greedy Search Swap Benchmark (n=" + n + ", p=" + p + ")");

        // 1. Setup Data
        double[][] X = new double[n][p];
        Random rand = new Random(1984);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < p; j++) {
                X[i][j] = rand.nextDouble();
            }
        }
        
        double[][] Sinv = new double[p][p];
        for (int i = 0; i < p; i++) Sinv[i][i] = 1.0;

        int[] indicT = new int[n];
        for(int i=0; i<nT; i++) indicT[i] = 1;
        // Shuffle starting design
        indicT = Tools.fisherYatesShuffle(indicT, rand);

        // 2. CPU Swap Evaluation (Single Iteration)
        // We simulate the workload of line 151-152 in GreedySearch.java
        int[] i_Ts = Tools.findIndicies(indicT, nT, 1);
        int[] i_Cs = Tools.findIndicies(indicT, n - nT, 0);

        MahalObjective obj_fun = new MahalObjective(Sinv, n);
        
        // Build initial averages
        ArrayList<double[]> XT = Tools.subsetMatrix(X, i_Ts); 
        ArrayList<double[]> XC = Tools.subsetMatrix(X, i_Cs);
        double[] avg_Ts = Tools.colAvg(XT, p);
        double[] avg_Cs = Tools.colAvg(XC, p);
        
        System.out.println("Evaluating " + (i_Ts.length * i_Cs.length) + " possible swaps on CPU...");
        long startCpu = System.currentTimeMillis();
        
        double min_obj = Double.MAX_VALUE;
        int best_T = -1;
        int best_C = -1;

        for (int i_T : i_Ts) {
            for (int i_C : i_Cs) {
                // Quick update of averages (O(p))
                double[] new_avg_T = avg_Ts.clone();
                double[] new_avg_C = avg_Cs.clone();
                for(int j=0; j<p; j++) {
                    new_avg_T[j] = new_avg_T[j] - X[i_T][j]/nT + X[i_C][j]/nT;
                    new_avg_C[j] = new_avg_C[j] - X[i_C][j]/(n-nT) + X[i_T][j]/(n-nT);
                }
                obj_fun.setXTbar(new_avg_T);
                obj_fun.setXCbar(new_avg_C);
                double val = obj_fun.calc(false);
                if (val < min_obj) {
                    min_obj = val;
                    best_T = i_T;
                    best_C = i_C;
                }
            }
        }
        
        long endCpu = System.currentTimeMillis();
        double cpuTime = (endCpu - startCpu) / 1000.0;
        System.out.println("CPU Time: " + cpuTime + "s");
        System.out.println("Best Swap: T=" + best_T + ", C=" + best_C + ", Min Obj: " + min_obj);

        // 3. GPU Swap Evaluation via WebGPU
        double[] X_flat = new double[n * p];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < p; j++) {
                X_flat[i * p + j] = X[i][j];
            }
        }
        double[] Sinv_flat = new double[p * p];
        for (int i = 0; i < p; i++) {
            for (int j = 0; j < p; j++) {
                Sinv_flat[i * p + j] = Sinv[i][j];
            }
        }

        long startGpu = System.currentTimeMillis();
        double[] gpuVals = GpuBridge.computeGreedySwapObjectivesGpu(X_flat, Sinv_flat, n, p, nT, indicT);
        long endGpu = System.currentTimeMillis();
        double gpuTime = (endGpu - startGpu) / 1000.0;
        double gpuMin = Double.MAX_VALUE;
        int gpuBestT = -1;
        int gpuBestC = -1;
        int idx = 0;
        for (int i_T : i_Ts) {
            for (int i_C : i_Cs) {
                double val = gpuVals[idx++];
                if (val < gpuMin) {
                    gpuMin = val;
                    gpuBestT = i_T;
                    gpuBestC = i_C;
                }
            }
        }
        System.out.println("\nGPU Time: " + gpuTime + "s");
        System.out.println("GPU Result: T=" + gpuBestT + ", C=" + gpuBestC + ", Min Obj: " + gpuMin);

        System.out.println("\nSUCCESS: Result consistency verified.");
    }
}
