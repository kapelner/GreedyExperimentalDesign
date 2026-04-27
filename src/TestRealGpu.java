import OptimalExperimentalDesign.GpuBridge;
import ObjectiveFunctions.KernelObjective;
import java.util.Random;

public class TestRealGpu {
    public static void main(String[] args) {
        int n = 5000;
        int r = 1000;
        int device_id = 1;

        System.out.println("Starting REAL GPU Benchmark from Java (n=" + n + ", r=" + r + ")");

        // 1. Setup Data
        Random rand = new Random(1984);
        double[][] Kgram = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = i; j < n; j++) {
                double val = rand.nextDouble();
                Kgram[i][j] = val; Kgram[j][i] = val;
            }
        }
        int[] indicT = new int[n];
        for (int i = 0; i < n; i++) indicT[i] = rand.nextBoolean() ? 1 : 0;

        // 2. CPU Calculation
        System.out.println("Running " + r + " evaluations on CPU...");
        long startCpu = System.currentTimeMillis();
        double cpuRes = 0;
        KernelObjective ko = new KernelObjective(Kgram);
        ko.setW(indicT);
        for(int i=0; i<r; i++) {
            ko.resetKernelSum();
            cpuRes = ko.calc(false);
        }
        long endCpu = System.currentTimeMillis();
        double cpuTime = (endCpu - startCpu)/1000.0;
        System.out.println("CPU Total Time: " + cpuTime + "s, Last Result: " + cpuRes);

        // 3. REAL GPU Calculation via JNI
        double[] W_flat_row_major = new double[r * n];
        for(int row=0; row<r; row++) {
            for(int i=0; i<n; i++) {
                W_flat_row_major[row * n + i] = (indicT[i] == 1 ? 1.0 : -1.0);
            }
        }
        double[] K_flat_row_major = new double[n * n];
        for(int i=0; i<n; i++) {
            for(int j=0; j<n; j++) {
                K_flat_row_major[i * n + j] = Kgram[i][j];
            }
        }

        GpuBridge gpu = new GpuBridge();
        System.out.println("Running " + r + " evaluations on GPU...");
        long startGpu = System.currentTimeMillis();
        double gpuRes = gpu.computeObjectiveGpu(W_flat_row_major, K_flat_row_major, r, n, device_id);
        long endGpu = System.currentTimeMillis();
        double gpuTime = (endGpu - startGpu)/1000.0;
        
        System.out.println("REAL GPU Total Time: " + gpuTime + "s, Result: " + gpuRes);
        
        if (Math.abs(cpuRes - gpuRes) < 1e-5) {
            System.out.println("SUCCESS: Results are consistent.");
            System.out.println("Real Speedup: " + (cpuTime / gpuTime) + "x");
        } else {
            System.out.println("FAILURE: Results differ! CPU=" + cpuRes + " GPU=" + gpuRes);
        }
    }
}
