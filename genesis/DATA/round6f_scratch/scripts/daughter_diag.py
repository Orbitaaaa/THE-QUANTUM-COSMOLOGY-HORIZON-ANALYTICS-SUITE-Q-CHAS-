import pickle, numpy as np, sys, glob, json
def diag(R, P, T, h, rN=np.sqrt(300.0), v=0.0015):
    ok = np.isfinite(R)
    fv = ok & (np.abs(P / v + 1) < 0.5); tv = ok & (np.abs(P / v - 1) < 0.5); big = ok & (R > 3 * rN)
    rows = []
    for i in range(R.shape[0]):
        js = np.where(fv[i] & big[i])[0]
        if not len(js):
            continue
        jl = js.max()
        after = np.where(ok[i] & (np.arange(R.shape[1]) > jl))[0]
        ex = after[tv[i, after] & big[i, after]] if len(after) else after
        jend = np.where(ok[i])[0].max()
        rows.append(dict(u=i * h, T_fv_end=float(T[i, jl]), r_at_fv_end=float(R[i, jl]), exit=bool(len(ex) > 0),
                         T_exit=float(T[i, ex[0]]) if len(ex) else None, r_max=float(np.nanmax(R[i, ok[i]])), r_end=float(R[i, jend]),
                         recollapse=bool(R[i, jend] < 0.05 * np.nanmax(R[i, ok[i]])), reaches_edge=bool(jend == R.shape[1] - 1)))
    return rows
if __name__ == "__main__":
    for f in sorted(glob.glob(sys.argv[1])):
        d = pickle.load(open(f, "rb"))
        rows = diag(d["R"], d["P"], d["T"], d["h"])
        ex = [r for r in rows if r["exit"]]
        rc = [r for r in rows if r["recollapse"]]
        print(f.split("/")[-1], "daughter rows", len(rows), "exiting", len(ex), "recollapsing", len(rc),
              "T_exit min/med", (min(r["T_exit"] for r in ex), float(np.median([r["T_exit"] for r in ex]))) if ex else None,
              "r_max of exiting rows", max(r["r_max"] for r in ex) if ex else None, "rows reaching edge", sum(r["reaches_edge"] for r in rows),
              "| T_exit_pred", d["out"]["T_exit_predicted"])
