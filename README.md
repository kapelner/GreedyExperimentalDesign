# GreedyExperimentalDesign

An R package for computing near optimal experimental designs for a two-arm experiment with covariates.

The repository contains a testscript as well for working through some examples.

This is joint work with Abba M. Krieger of the Wharton School of the University of Pennsylvania and David Azriel of The Technion.

To load the package, make sure `rJava' is installed and properly configured! There are some issues with the latest Java 10. Make sure you read about this online. Then:

	options(java.parameters = "-Xmx4000m")
	library(GreedyExperimentalDesign)
	
And if you want to use the optimization feature via `Gurobi', install Gurobi first and then run something like:

	.jaddClassPath("/gurobi752/win64/lib/gurobi.jar")
