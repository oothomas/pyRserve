# --- NumPy 2.0 compat shim for pyRserve (must come before importing pyRserve) ---
import numpy as _np
from types import SimpleNamespace as _SimpleNamespace
if not hasattr(_np, "string_"):  _np.string_  = _np.bytes_
if not hasattr(_np, "unicode_"): _np.unicode_ = str
if not hasattr(_np, "int"):      _np.int      = int
if not hasattr(_np, "bool"):     _np.bool     = bool
if not hasattr(_np, "float"):    _np.float    = float
if not hasattr(_np, "object"):   _np.object   = object
if not hasattr(_np, "compat"):   _np.compat   = _SimpleNamespace(long=int)
elif not hasattr(_np.compat, "long"):
    _np.compat.long = int
# -------------------------------------------------------------------------------

import pyRserve
import pandas as pd
import numpy as np
import slicer
import Support.gpa_lib as gpa_lib

conn = pyRserve.connect()
conn.r.ls()

# -------------------------------------------------------------------------
# Replace CSV I/O with direct read from SlicerMorph GPA widget
# -------------------------------------------------------------------------
gpaWidget = slicer.modules.gpa.widgetRepresentation().self()
lmData = gpaWidget.LM  # SlicerMorph's LMData object (post-GPA data)

# Aligned landmarks are in lmData.lm with shape (p, 3, n)
arr = np.asarray(lmData.lm)  # p × 3 × n
assert arr.ndim == 3 and arr.shape[1] == 3, f"Unexpected LM.lm shape: {arr.shape}"
p, _, n = arr.shape

# Build coords as n × (3p) in x1,y1,z1,x2,y2,z2,... order
coords_mat = arr.transpose(2, 0, 1).reshape(n, 3 * p, order="C")

# Centroid sizes (vector of length n). Module spells it 'centriodSize'.
size = np.asarray(lmData.centriodSize, dtype=float).reshape(-1)
if size.shape[0] != n:
    raise RuntimeError(f"centroid size length {size.shape[0]} != number of specimens {n}")

# Send to R as before
conn.r.coords = coords_mat
conn.r.size = size

conn.eval('require(geomorph)')

# arrayspecs expects coords as n × (k*p); use the real p from lmData
conn.eval(f'arr=arrayspecs(coords,p={p},k=3)')

#only variable in this linear model is size against shape (i.e., LM coordinates)
#in reality these will be provided by the user inside the GPA module
#things such as Sex, Genotype, locomotion, ecology etc

conn.voidEval('gdf=geomorph.data.frame(Size=size,Coords=arr)')

#here is the linear model 
#in the actual GPA module, this will be a user input based on the covariates they created
#size is always an option. It is derived from LM coordinates. 
conn.voidEval('mod=as.formula(Coords~Size)')

#other example models
# conn.voidEval('mod=as.formula(Coords~Size*Sex+Genotype))

#that in R returns a two matrix, first row is the intercept, and the second is the coefs of size
#each subsequent variable in the formula will be other rows.
conn.voidEval('outlm=procD.lm(mod, data=gdf)')
model=conn.eval('outlm')
#conn.eval('remove(list=ls())')
#conn.shutdown()
#to shutdown Rserve. Otherwise it will stay open and contain all variable from the current session).


# --- Extract the regression coefficients from the R model ---
# modelCoeffs is an AttrArray of shape (2, 3p) with row 0 = intercept, row 1 = size slope
modelCoeffs = model[2]
intercept_row = modelCoeffs[0]  # shape (3p,)
slope_row     = modelCoeffs[1]  # shape (3p,)

print("modelCoeffs.shape", modelCoeffs.shape)
print("intercept_row.shape", intercept_row.shape)
print("slope_row.shape", slope_row.shape)

print("Before adding slope:")
print("  pcNumber:", gpaWidget.pcNumber)
print("  lmData.vec:", lmData.vec.shape)
print("  lmData.val:", lmData.val.shape)

# --- Verify number of landmarks (p) from lmData (already computed above) ---
print("GPA says #landmarks p =", p, "; slope_row =", slope_row.shape)

# --- Reshape the slope row into (p, 3) and then flatten to (3p,) ---
# Use a high multiplier to exaggerate the effect; adjust the value if needed.
size_slope_flat = slope_row.reshape(p, 3, order='C').ravel(order='C')
size_slope_flat *= 1000  # Increase multiplier

# --- Append the new PC vector ---
lmData.vec = np.column_stack((lmData.vec, size_slope_flat))
dummy_eigenvalue = np.median(lmData.val)
lmData.val = np.append(lmData.val, dummy_eigenvalue)

print("After adding slope:")
print("  lmData.vec:", lmData.vec.shape)
print("  lmData.val:", lmData.val.shape)

# --- Sync pcNumber and update the GPA widget ---
gpaWidget.pcNumber = lmData.vec.shape[1]

lmData.sortedEig = gpa_lib.pairEig(lmData.val, lmData.vec)
gpaWidget.updateList()
print("All done. New PC has been added.")

# --- Debug: print the norm of the new PC vector vs PC1 ---
newPC = lmData.vec[:, -1]
pc1 = lmData.vec[:, 0]
print("Norm of new PC (size slope):", np.linalg.norm(newPC))
print("Norm of PC1:", np.linalg.norm(pc1))
