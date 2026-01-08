# GreedyExperimentalDesign

An R package for computing near optimal experimental designs for a two-arm experiment with covariates. This is joint work with Abba M. Krieger of the Wharton School of the University of Pennsylvania and David Azriel of The Technion.

To load the package, make sure `rJava` is installed and properly configured! Then:

	options(java.parameters = "-Xmx5g") #or however much memory you wish to allocate
	library(GreedyExperimentalDesign)
	
If you want to use the numerical optimization design, download [https://www.gurobi.com/downloads/gurobi-software/](Gurobi) and register your machine with your license via `grbgetkey`. Then, this download comes with its own proprietary R package and you need to do a local install via `R CMD INSTALL`.

## Manual cleanup before install

If `R CMD INSTALL GreedyExperimentalDesign` fails with architecture mismatch errors (for example, "not a valid Win32 application"), remove stale compiled artifacts in the package `src/` directory and reinstall.

On macOS/Linux (from `GreedyExperimentalDesign/`):

```
rm -f src/*.o src/*.so src/*.dll src/*.dylib
```

On Windows (from `GreedyExperimentalDesign\`):

```
del /Q src\\*.o src\\*.so src\\*.dll src\\*.dylib
```

## Mahalanobis benchmark

The following benchmark demonstrates this package's many designs. It draws a fixed covariate matrix and plots the distribution of Mahalanobis distances for the candidate designs each method returns. Lower proportional Mahalanobis distance means better balance, so methods whose histograms are shifted further left tend to produce better-balanced designs on this dataset. The vertical lines show each method's average proportional Mahalanobis distance. Complete randomization provides a baseline; methods that concentrate mass well to the left of that baseline are producing stronger balance improvements under the same objective.

![mahal_dist_all_searches.png 1080x720](mahal_dist_all_searches.png)

To generate the figure above, run:

	Rscript GreedyExperimentalDesign/benchmarks/plot_obj_value_comparison_all_searches.R
