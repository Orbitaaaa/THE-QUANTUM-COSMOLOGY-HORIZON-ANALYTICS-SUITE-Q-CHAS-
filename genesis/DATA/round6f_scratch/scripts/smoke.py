import sys, time, json, types
sys.argv = ["x", "--figure-only"]
src = open("round6f.py").read()
src = src.split("# ------------------------------------------------------------------ run")[0]
ns = {"__name__": "smoke"}
exec(compile(src, "round6f_head", "exec"), ns)
t = time.time(); print(ns["threshold_job"]((4.0, 0.3, 0.999, 30.0)), round(time.time() - t, 1), flush=True)
t = time.time(); print(ns["continuous_job"]((2.0, 1e-6, 1.5e-9)), round(time.time() - t, 1), flush=True)
