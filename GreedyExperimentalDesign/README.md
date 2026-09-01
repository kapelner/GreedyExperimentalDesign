# GreedyExperimentalDesign

[![R-universe version](https://kapelner.r-universe.dev/GreedyExperimentalDesign/badges/version)](https://kapelner.r-universe.dev/GreedyExperimentalDesign)
[![R-universe checks](https://kapelner.r-universe.dev/GreedyExperimentalDesign/badges/checks)](https://kapelner.r-universe.dev/GreedyExperimentalDesign)
[![R-CMD-check](https://github.com/kapelner/GreedyExperimentalDesign/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/kapelner/GreedyExperimentalDesign/actions/workflows/R-CMD-check.yaml)

`GreedyExperimentalDesign` computes covariate-balanced treatment allocations for two-arm experiments using greedy optimization, rerandomization, matching, Karp's method, exhaustive search, Hadamard designs, and multiple-kernel objectives.

## Installation

Install the current development version from Adam Kapelner's R-universe:

```r
install.packages(
  "GreedyExperimentalDesign",
  repos = c(
    kapelner = "https://kapelner.r-universe.dev",
    CRAN = "https://cloud.r-project.org"
  )
)
```

The package requires Java 21 or newer. Set JVM options before loading the package. The preview flag is needed only when using the optional Java WebGPU backend, but is safe to set for CPU-only use as well.

```r
options(java.parameters = c(
  "-Xmx10g",
  "--enable-native-access=ALL-UNNAMED",
  "--enable-preview"
))
library(GreedyExperimentalDesign)
```

See the [repository README](https://github.com/kapelner/GreedyExperimentalDesign#readme) for GPU setup and extended usage instructions.
