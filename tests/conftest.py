import os

# Dense 957x957 linear algebra is fastest single-threaded on shared CI runners;
# multi-threaded OpenBLAS was measured ~60x slower under contention.
for _var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_var, "1")
