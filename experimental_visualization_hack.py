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

import os
import pyRserve
import pandas as pd
import numpy as np
import slicer
import Support.gpa_lib as gpa_lib

# -------------------- USER TUNABLES --------------------
RSERVE_HOST   = '127.0.0.1'
RSERVE_PORT   = 6311
RSERVE_TO     = 3            # unused now; kept for consistency
MAGNIFY       = 1000.0       # exaggeration for visibility of added vectors
TARGET_LEVEL  = None         # e.g., "Pan" to force which covariate level to use; None = first non-reference
# -------------------------------------------------------

def _connect_rserve(host=RSERVE_HOST, port=RSERVE_PORT, timeout=RSERVE_TO):
    try:
        return pyRserve.connect(host='127.0.0.1', port=6311)
    except (ConnectionRefusedError, OSError, Exception) as e:
        raise RuntimeError(f"Cannot reach Rserve at 127.0.0.1:6311: {e}")

def _find_covariates_path(gpaWidget):
    # Prefer explicit selection in the GPA UI
    try_paths = []
    try:
        p = str(gpaWidget.selectCovariatesText.text)
        if p:
            try_paths.append(p)
    except Exception:
        pass
    # Fall back to the standard output copy (if GPA wrote it out)
    try:
        p2 = os.path.join(gpaWidget.outputFolder, "covariateTable.csv")
        try_paths.append(p2)
    except Exception:
        pass
    # Return the first existing
    for p in try_paths:
        if p and os.path.isfile(p):
            return p
    return None

def _load_covariate_series(cov_path, files):
    """
    Read covariate CSV, align to specimen order in `files`.
    Expects column 0 = ID, column 1 = first covariate (categorical).
    Returns (cov_name, cov_series_aligned, unique_levels)
    """
    df = pd.read_csv(cov_path)
    if df.shape[1] < 2:
        raise ValueError(f"Covariate table {cov_path} must have at least 2 columns (ID + one covariate).")
    id_col = df.columns[0]
    cov_name = df.columns[1]
    # Index by ID and reindex to files order
    df = df.set_index(id_col)
    missing = [f for f in files if f not in df.index]
    if missing:
        raise ValueError(f"IDs missing in covariate table (first 10 shown): {missing[:10]}")
    cov = df[cov_name].astype(str).reindex(files)
    levels = list(pd.Categorical(cov).categories)
    if len(levels) < 2:
        raise ValueError(f"Covariate '{cov_name}' must have >=2 levels; found {levels}")
    return cov_name, cov.tolist(), levels

def _rowpick(names, key):
    # Return index of the row whose name equals `key`; fallback to contains
    for i, nm in enumerate(names):
        if nm == key:
            return i
    for i, nm in enumerate(names):
        if key in nm:
            return i
    return None

# Launch R and, run the following:
# library(Rserve); Rserve(debug=TRUE, args="--no-save")

# -------------------- MAIN --------------------
# 1) Hook into GPA widget and grab aligned landmark data
gpaWidget = slicer.modules.gpa.widgetRepresentation().self()
lmData = gpaWidget.LM  # SlicerMorph's LMData (post-GPA)
arr = np.asarray(lmData.lm)  # shape (p, 3, n)
assert arr.ndim == 3 and arr.shape[1] == 3, f"Unexpected LM.lm shape: {arr.shape}"
p, _, n = arr.shape

# Build coords as n × (3p): x1,y1,z1, x2,y2,z2, …
coords_mat = arr.transpose(2, 0, 1).reshape(n, 3 * p, order="C")

# Centroid sizes (vector of length n) — attribute name is 'centriodSize'
size_vec = np.asarray(lmData.centriodSize, dtype=float).reshape(-1)
if size_vec.shape[0] != n:
    raise RuntimeError(f"centroid size length {size_vec.shape[0]} != number of specimens {n}")

# Specimen ID order (must match coords/size order)
if hasattr(gpaWidget, "files") and isinstance(gpaWidget.files, (list, tuple)) and len(gpaWidget.files) == n:
    files = list(gpaWidget.files)
else:
    files = [f"spec_{i+1}" for i in range(n)]

# 2) Optional covariates
covariates_path = _find_covariates_path(gpaWidget)
have_cov = False
cov_name = None
cov_list = None
cov_levels = None
if covariates_path:
    try:
        cov_name, cov_list, cov_levels = _load_covariate_series(covariates_path, files)
        have_cov = True
        print(f"[Covariates] Using '{cov_name}' from: {covariates_path}")
        print(f"[Covariates] Levels: {cov_levels}")
    except Exception as e:
        print(f"[Covariates] Skipping covariate model ({e}). Proceeding with Size-only.")

# 3) Connect to Rserve and fit model(s)
conn = _connect_rserve()

# Ship data to R
conn.r.coords = coords_mat
conn.r.size   = size_vec

# Load needed pkgs and normalize inputs on the R side
conn.eval('suppressPackageStartupMessages({'
          'if (!requireNamespace("geomorph", quietly=TRUE)) stop("Package geomorph not installed");'
          'if (!requireNamespace("RRPP", quietly=TRUE)) stop("Package RRPP not installed");'
          'library(geomorph); library(RRPP)'
          '})')

# Ensure coords is a numeric matrix before arrayspecs()
conn.eval('coords <- as.matrix(coords); storage.mode(coords) <- "double"')
conn.eval(f'arr <- arrayspecs(coords, p={p}, k=3)')

if have_cov:
    # Provide covariate as Python list -> coerce to character vector -> factor
    conn.r.cov = cov_list
    conn.eval('cov <- factor(as.character(unlist(cov)))')
    conn.eval('gdf <- geomorph.data.frame(Size=as.numeric(size), Covariate=cov, Coords=arr)')
    conn.eval('mod <- as.formula("Coords ~ Size * Covariate")')
    conn.eval('outlm <- procD.lm(mod, data=gdf)')
    # Extract coefficients with names and levels explicitly
    conn.eval('coef_mat <- outlm$coefficients')
    conn.eval('coef_names <- rownames(coef_mat)')
    conn.eval('levs <- levels(gdf$Covariate)')
    coef_mat   = np.asarray(conn.eval('coef_mat'))         # shape: (#terms, 3p)
    coef_names = list(conn.eval('as.character(coef_names)'))
    cov_levels_r = list(conn.eval('as.character(levs)'))
else:
    conn.eval('gdf <- geomorph.data.frame(Size=as.numeric(size), Coords=arr)')
    conn.eval('mod <- as.formula("Coords ~ Size")')
    conn.eval('outlm <- procD.lm(mod, data=gdf)')
    conn.eval('coef_mat <- outlm$coefficients')
    conn.eval('coef_names <- rownames(coef_mat)')
    coef_mat   = np.asarray(conn.eval('coef_mat'))
    coef_names = list(conn.eval('as.character(coef_names)'))

# We’re done with R for now
try:
    conn.shutdown()
except Exception:
    pass

# 4) Pull out slopes and append as vectors to GPA
def _flat_3p(row_3p):
    """row_3p: (3p,) -> flattened in (x1,y1,z1,...), already OK. Provide as (3p,) contiguous."""
    v = np.asarray(row_3p).reshape(p, 3, order='C').ravel(order='C')
    return v

new_vectors = []
new_names   = []

# Always include Size (allometry) slope
idx_size = _rowpick(coef_names, "Size")
if idx_size is None:
    raise RuntimeError(f"Could not find 'Size' in coefficient names: {coef_names}")
size_vec_3p = _flat_3p(coef_mat[idx_size, :]) * MAGNIFY
new_vectors.append(size_vec_3p)
new_names.append("Allometry (Size)")

if have_cov:
    # Decide which level to use (first non-reference by default or user-specified)
    ref = cov_levels_r[0]
    candidates = cov_levels_r[1:] if len(cov_levels_r) > 1 else []
    if TARGET_LEVEL is not None:
        if TARGET_LEVEL not in cov_levels_r:
            raise RuntimeError(f"TARGET_LEVEL '{TARGET_LEVEL}' not in levels {cov_levels_r}")
        if TARGET_LEVEL == ref:
            alts = [lv for lv in cov_levels_r if lv != ref]
            if not alts:
                raise RuntimeError(f"No alternative level to contrast with reference '{ref}'.")
            chosen = alts[0]
        else:
            chosen = TARGET_LEVEL
    else:
        chosen = candidates[0] if candidates else ref

    # Main effect row name looks like Covariate<Level> (treatment contrasts)
    main_key = f"Covariate{chosen}"
    idx_main = _rowpick(coef_names, main_key) or _rowpick(coef_names, f"Covariate={chosen}")
    if idx_main is None:
        raise RuntimeError(f"Could not find main-effect row for level '{chosen}' in {coef_names}")

    cov_main_vec = _flat_3p(coef_mat[idx_main, :]) * MAGNIFY
    new_vectors.append(cov_main_vec)
    new_names.append(f"{cov_name}: {chosen} vs {ref}")

    # Interaction row name looks like Size:Covariate<Level>
    inter_key = f"Size:Covariate{chosen}"
    idx_inter = _rowpick(coef_names, inter_key) or _rowpick(coef_names, f"Size:Covariate={chosen}")
    if idx_inter is None:
        raise RuntimeError(f"Could not find interaction row for '{chosen}' in {coef_names}")

    inter_vec = _flat_3p(coef_mat[idx_inter, :]) * MAGNIFY
    new_vectors.append(inter_vec)
    new_names.append(f"Size×{cov_name}: {chosen}")

# 5) Append to LMData and refresh UI
for v in new_vectors:
    lmData.vec = np.column_stack((lmData.vec, v))
    lmData.val = np.append(lmData.val, 1.0)  # dummy eigenvalue for bookkeeping

# Sync & re-pair
gpaWidget.pcNumber = lmData.vec.shape[1]
lmData.sortedEig = gpa_lib.pairEig(lmData.val, lmData.vec)
gpaWidget.updateList()

# Recompute scatter data for ALL components/vectors
numPC  = gpaWidget.pcNumber
n_spec = lmData.lm.shape[2]
gpaWidget.scatterDataAll = np.zeros((n_spec, numPC))
for i in range(numPC):
    scores = gpa_lib.plotTanProj(lmData.lm, lmData.sortedEig, i, 1)  # (n × ?), scores in col 0
    gpaWidget.scatterDataAll[:, i] = scores[:, 0]

# Rename the last k entries (k = len(new_names)) everywhere so they’re not mislabeled as PCs
k = len(new_names)
if k > 0:
    try:
        # Internal list used by UI
        if hasattr(gpaWidget, "PCList") and len(gpaWidget.PCList) >= k:
            for i, nm in enumerate(new_names[::-1]):
                gpaWidget.PCList[-(i+1)] = nm

        # Rename in combo boxes
        total = gpaWidget.XcomboBox.count
        if total >= k:
            for i, nm in enumerate(new_names[::-1]):
                idx = total - 1 - i
                for cb in [gpaWidget.XcomboBox, gpaWidget.YcomboBox,
                           gpaWidget.vectorOne, gpaWidget.vectorTwo, gpaWidget.vectorThree]:
                    if cb.count > idx:
                        cb.setItemText(idx, nm)

        # Refresh slider-group dropdown
        if hasattr(gpaWidget, "slider1") and hasattr(gpaWidget, "PCList"):
            gpaWidget.slider1.populateComboBox(gpaWidget.PCList)
            if gpaWidget.slider1.comboBox.count >= k:
                for i, nm in enumerate(new_names[::-1]):
                    idx = gpaWidget.slider1.comboBox.count - 1 - i
                    gpaWidget.slider1.comboBox.setItemText(idx, nm)
    except Exception as e:
        print("Rename UI entries warning:", e)

print(f"Added {len(new_names)} vector(s): {', '.join(new_names)}")
print("Done.")