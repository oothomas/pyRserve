# pyRserve
Minimal pyRserve example

## Prerequsites

* Install Rserve and geomorph packages to your local R installation via `install.packages(c('Rserve', 'geomorph'))`
* Install pyRserve and pandas packages to your Slicer via `pip_install("pyRserve pandas")`

## Steps

First, to start Rserve session (Create a separate standalone R via commandline, not Rstudio) type:
```
library(geomorph)
library(Rserve)
run.Rserve()
```

Updated the hard coded path in `size_regression.py` to where you have cloned this repo:
https://github.com/johnbradley/pyRserve/blob/f3d85b47cf2acc9eae6ac7906f405ea1097c28c9/size_regression.py#L9

Then, copy and paste the size_regression.py into Python console of Slicer.

You can compare the contents of the summary variable in python to the same R code executed in Rstudio. 
