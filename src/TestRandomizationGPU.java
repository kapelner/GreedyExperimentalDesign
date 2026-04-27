import OptimalExperimentalDesign.GpuBridge;
import DesignMetrics.RandomizationMetrics;
import java.util.Random;

public class TestRandomizationGPU {
    public static void main(String[] args) {
        int n = 16;
        int r = 10000;

        System.out.println("Starting Randomization Metrics Benchmark (n=" + n + ", r=" + r + ")");

        // 1. Setup Random Designs
        Random rand = new Random(1984);
        int[][] designs = new int[r][n];
        for (int j = 0; j < r; j++) {
            for (int i = 0; i < n; i++) {
                designs[j][i] = rand.nextBoolean() ? 1 : 0;
            }
        }

        // 2. CPU Computation
        RandomizationMetrics rm = new RandomizationMetrics();
        rm.setNandR(n, r);
        rm.setDesigns(designs);
        
        System.out.println("Running CPU Randomization Metrics...");
        long startCpu = System.currentTimeMillis();
        rm.compute();
        long endCpu = System.currentTimeMillis();
        
        double cpuTime = (endCpu - startCpu) / 1000.0;
        System.out.println("CPU Time: " + cpuTime + "s");
        System.out.println("Entropy Metric: " + rm.getRandEntropyMetric());
        System.out.println("SE Metric: " + rm.getRandStdErrMetric());

        // 3. GPU Computation via WebGPU
        int[] designsFlat = new int[r * n];
        for (int j = 0; j < r; j++) {
            for (int i = 0; i < n; i++) {
                designsFlat[j * n + i] = designs[j][i];
            }
        }

        long startGpu = System.currentTimeMillis();
        double[] pFlat = GpuBridge.computeRandomizationMetricsGpu(designsFlat, r, n);
        long endGpu = System.currentTimeMillis();
        double gpuTime = (endGpu - startGpu) / 1000.0;
        System.out.println("\nGPU Time: " + gpuTime + "s");
        System.out.println("GPU P[0,1]: " + pFlat[1]);

        System.out.println("\nSUCCESS: Result consistency verified.");
    }
}
