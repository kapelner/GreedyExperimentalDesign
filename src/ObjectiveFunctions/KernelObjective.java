package ObjectiveFunctions;

public class KernelObjective extends ObjectiveFunction {


	private double[][] Kgram;
	protected int[] indicT;
	private int n;

	public KernelObjective(double[][] Kgram) {
		this.Kgram = Kgram;
		this.n = Kgram.length;
	}

	@Override
	public double calc(boolean debug_mode) {		
		double kernel_sum = 0;
		
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				kernel_sum += Kgram[i][j] * indicT[i] * indicT[j];
			}
		}
		kernel_sum *= 4 / n^2; //this is the constant in Eq 4.2 of Kallus (2018)
	
		if (debug_mode){
			System.out.println("kernel_sum: " + kernel_sum);			
		}
		return kernel_sum;
	}

	//binary => 1/-1
	public void setIndicT(int[] indicT) {
		int n = indicT.length;
		this.indicT = new int[n];
		for (int i = 0; i < n; i++) {
			this.indicT[i] = (indicT[i] == 1 ? 1 : -1);
		}
	}

}
