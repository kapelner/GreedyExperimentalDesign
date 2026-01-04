# GreedyExperimentalDesign

An R package for computing near optimal experimental designs for a two-arm experiment with covariates.

The repository contains a testscript as well for working through some examples.

This is joint work with Abba M. Krieger of the Wharton School of the University of Pennsylvania and David Azriel of The Technion.

## Mahalanobis benchmark

![mahal_dist_all_searches.png 1080x720](mahal_dist_all_searches.png)

To generate the figure above, run:

	Rscript GreedyExperimentalDesign/benchmarks/plot_obj_value_comparison_all_searches.R

This benchmark draws a fixed covariate matrix, runs several search methods, and plots the distribution of Mahalanobis distances for the candidate designs each method returns. Lower Mahalanobis distance means better balance, so methods whose histograms are shifted further left and whose mean lines are further left tend to produce better-balanced designs on this dataset. The vertical lines show each method's mean objective value. Complete randomization provides a baseline; methods that concentrate mass well to the left of that baseline are producing stronger balance improvements under the same objective.

To load the package, make sure `rJava' is installed and properly configured! There are some issues with the latest Java 10. Make sure you read about this online. Then:

	options(java.parameters = "-Xmx4000m")
	library(GreedyExperimentalDesign)
	
And if you want to use the optimization feature via `Gurobi', install Gurobi first and then run something like:

	.jaddClassPath("/gurobi801/win64/lib/gurobi.jar")
	
Note that you have to register your machine and run `grbgetkey` with the proper license number and ensure the `.lic` file is installed. Expect `ClassNotFoundException`s if this is not setup properly!
