"""BenchFEA: 3D frame analysis engine and structural-modelling benchmark."""
import os as _os

# The matrices of this package are small and dense; a multi-threaded BLAS only costs time on them
# (measured: orders of magnitude on some machines).  The BLAS libraries read these variables when
# they are loaded, so they are set here, before anything of the package imports numpy; a value set
# by the caller is kept, and a process that imported numpy before benchfea is not affected.
for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    _os.environ.setdefault(_var, "1")
del _os, _var

__version__ = "0.1.0"
