# pyRserve
Minimal pyRserve example

## Prerequsites

* Install Rserve package to your local R installation
* Install pyRserve package to your Slicer via `pip_install("pyRserve")`


First, start R session (not Rstudio, the commandline R. It will block the UI and the terminal) via
```
library(Rserve)
run.Rserve()
```

Then, copy and paste the size_regression.py into Python console of Slicer.

You can compare the contents of the summary variable in python to the same R code executed in Rstudio. 
