# pyRserve
Minimal pyRserve example

## Prerequsites

* Install Rserve and geomorph packages to your local R installation. 
* Install pyRserve package to your Slicer via `pip_install("pyRserve")`


First, to start Rserve session (Create a separate standalone R via commandline, not Rstudio) type:
```
library(Rserve)
run.Rserve()
```

Then, copy and paste the size_regression.py into Python console of Slicer.

You can compare the contents of the summary variable in python to the same R code executed in Rstudio. 
