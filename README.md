# SpiceFP

**Sparse and Structured Procedure to Identify Combined Effects of Functional Predictors**

A set of functions allowing to implement the **spiceFP** approach which is 
iterative. It involves transformation of functional predictors into several 
candidate explanatory matrices (based on contingency tables), to which 
edge matrices with contiguity constraints are associated.

**Generalized Fused Lasso regression** are performed in order to identify the 
best candidate matrix, the best class intervals and related coefficients at 
each iteration. 

The approach is stopped when the maximal number of iterations is reached or 
when retained coefficients are zeros. Supplementary functions allow to get 
coefficients of any candidate matrix or mean of coefficients of many candidates.

# Installation

To install the **SpiceFP** package, the easiest is to install it directly 
from GitHub. Open an R session and run the following commands:

```R
library(remotes) 
install_github("giraultg/SpiceFP", build_vignettes=TRUE)
```

# Usage

Once the package is installed on your computer, it can be loaded into a R session:

```R
library(spiceFP)
help(package="SpiceFP")
```

# Citation

As a lot of time and effort were spent in creating the **SpiceFP** method, 
please cite it when using it for data analysis:

https://hal.archives-ouvertes.fr/hal-03298977

You should also cite the **SpiceFP** package:

```R
citation("spiceFP")
```

See also citation() for citing R itself.

# References

1. Taylor B. Arnold and Ryan J. Tibshirani (2020). genlasso: Path Algorithm for
  Generalized Lasso Problems. R package version 1.5.
  https://CRAN.R-project.org/package=genlasso
  
