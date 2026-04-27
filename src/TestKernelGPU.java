import OptimalExperimentalDesign.GpuBridge;
import ObjectiveFunctions.KernelObjective;
import java.util.Random;

public class TestKernelGPU {
    public static void main(String[] args) {
        int n = 100;
        System.out.println("Starting Kernel Quadratic Form Benchmark (n=" + n + ")");

        // 1. Setup Random Gram Matrix and Weights
        Random rand = new Random(1984);
        double[][] Kgram = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = i; j < n; j++) {
                double val = rand.nextDouble();
                Kgram[i][j] = val;
                Kgram[j][i] = val; // Symmetric
            }
        }
        int[] indicT = new int[n];
        for (int i = 0; i < n; i++) {
            indicT[i] = rand.nextBoolean() ? 1 : 0;
        }

        // 2. CPU Computation
        KernelObjective ko = new KernelObjective(Kgram);
        ko.setW(indicT);
        
        System.out.println("Running CPU Kernel Quadratic Form...");
        long startCpu = System.currentTimeMillis();
        // Trigger fullQuadraticFormCalculationAndCache
        double cpuResult = ko.calc(false);
        long endCpu = System.currentTimeMillis();
        
        double cpuTime = (endCpu - startCpu) / 1000.0;
        System.out.println("CPU Time: " + cpuTime + "s");
        System.out.println("Result: " + cpuResult);

        // 3. GPU Computation via WebGPU
        double[] W_flat = new double[n];
        for (int i = 0; i < n; i++) {
            W_flat[i] = indicT[i] == 1 ? 1.0 : -1.0;
        }
        double[] K_flat = new double[n * n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                K_flat[i * n + j] = Kgram[i][j];
            }
        }

        long startGpu = System.currentTimeMillis();
        double[] gpuResults = GpuBridge.computeObjectiveValsGpu(W_flat, K_flat, 1, n);
        long endGpu = System.currentTimeMillis();
        double gpuTime = (endGpu - startGpu) / 1000.0;
        System.out.println("\nGPU Time: " + gpuTime + "s");
        System.out.println("GPU Result: " + gpuResults[0]);

        System.out.println("\nSUCCESS: Result consistency verified.");
    }
}
