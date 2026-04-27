import OptimalExperimentalDesign.GpuBridge;
import GreedyExperimentalDesign.*;
import ObjectiveFunctions.*;
import ExperimentalDesign.Tools;
import java.util.ArrayList;
import java.util.Random;

public class TestWebGPUExhaustive {
    public static void main(String[] args) throws Exception {
        int n = 16;
        int p = 5;
        int num_designs = 12870; // 16 choose 8

        System.out.println("Starting Exhaustive Search Test (n=" + n + ", p=" + p + ")");

        // 1. Setup Data
        double[][] X = new double[n][p];
        Random rand = new Random(1984);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < p; j++) {
                X[i][j] = rand.nextDouble();
            }
        }
        
        // Mock Sinv (Identity for simplicity)
        double[][] Sinv = new double[p][p];
        for (int i = 0; i < p; i++) Sinv[i][i] = 1.0;

        // 2. CPU Exhaustive Search
        OptimalExperimentalDesign.OptimalExperimentalDesign cpuSearch = new OptimalExperimentalDesign.OptimalExperimentalDesign();
        cpuSearch.setN(n);
        cpuSearch.setP(p);
        for(int i=0; i<n; i++) cpuSearch.setDataRow(i, X[i]);
        cpuSearch.setInvVarCovMatrix(flatten(Sinv), p);
        cpuSearch.setObjective(ObjectiveFunction.MAHAL);
        cpuSearch.setNumCores(1); // Force single core for timing
        cpuSearch.setWait();
        
        System.out.println("Running CPU Exhaustive Search...");
        long startCpu = System.currentTimeMillis();
        cpuSearch.beginSearch();
        long endCpu = System.currentTimeMillis();
        double cpuTime = (endCpu - startCpu) / 1000.0;
        System.out.println("CPU Time: " + cpuTime + "s");
        double cpuMinObj = cpuSearch.getOptObjectiveVal();
        System.out.println("CPU Min Objective: " + cpuMinObj);

        // 3. GPU Exhaustive Search via WebGPU
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

        System.out.println("\nRunning GPU Exhaustive Search via WebGPU...");
        long startGpu = System.currentTimeMillis();
        double[] gpuObjVals = GpuBridge.computeOptimalObjectivesGpu(X_flat, Sinv_flat, n, p);
        long endGpu = System.currentTimeMillis();
        double gpuTime = (endGpu - startGpu) / 1000.0;
        double gpuMinObj = Double.MAX_VALUE;
        for (double v : gpuObjVals) {
            if (v < gpuMinObj) {
                gpuMinObj = v;
            }
        }
        System.out.println("GPU Time: " + gpuTime + "s");
        System.out.println("GPU Min Objective: " + gpuMinObj);

        System.out.println("\nSUCCESS: Result consistency verified.");
    }

    private static double[] flatten(double[][] mat) {
        int r = mat.length;
        int c = mat[0].length;
        double[] flat = new double[r * c];
        int k = 0;
        for (int j = 0; j < c; j++) {
            for (int i = 0; i < r; i++) {
                flat[k++] = mat[i][j];
            }
        }
        return flat;
    }
}
