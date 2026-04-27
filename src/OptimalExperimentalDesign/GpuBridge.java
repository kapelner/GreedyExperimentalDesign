package OptimalExperimentalDesign;

public class GpuBridge {
    public static double[] computeObjectiveValsGpu(double[] W_flat, double[] K_flat, int r, int n) {
        return WebGpuPanama.computeObjectiveVals(W_flat, K_flat, r, n);
    }

    public static double[] computeRandomizationMetricsGpu(int[] designs_flat, int r, int n) {
        return WebGpuPanama.computeRandomizationMetrics(designs_flat, r, n);
    }

    public static double[] computeGreedySwapObjectivesGpu(
        double[] X_flat,
        double[] Sinv_flat,
        int n,
        int p,
        int nT,
        int[] indicT
    ) {
        return WebGpuPanama.computeGreedySwapObjectives(X_flat, Sinv_flat, n, p, nT, indicT);
    }

    public static double[] computeOptimalObjectivesGpu(double[] X_flat, double[] Sinv_flat, int n, int p) {
        return WebGpuPanama.computeOptimalObjectives(X_flat, Sinv_flat, n, p);
    }

    public double computeObjectiveGpu(double[] W_flat, double[] K_flat, int r, int n, int device) {
        double[] vals = computeObjectiveValsGpu(W_flat, K_flat, r, n);
        return vals.length == 0 ? Double.NaN : vals[0];
    }

    public double computeObjectiveGpu(double[] W_flat, double[] K_flat, int r, int n) {
        return computeObjectiveGpu(W_flat, K_flat, r, n, 0);
    }
}
