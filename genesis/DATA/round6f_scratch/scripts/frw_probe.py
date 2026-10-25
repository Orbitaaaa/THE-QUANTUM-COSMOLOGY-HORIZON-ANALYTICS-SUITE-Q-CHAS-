import numpy as np
S = "/tmp/claude-0/-home-user-THE-QUANTUM-COSMOLOGY-HORIZON-ANALYTICS-SUITE-Q-CHAS-/5d1f4f78-e2f6-5e44-a419-ab5b29ff8104/scratchpad/r7"
d = np.load(S + "/route_full_condensate.npz")
R, P, T, n = d["r"].astype(float), d["phi"], d["T"].astype(float), d["n"]
v, rN = 0.0015, np.sqrt(300.0)
ok = np.isfinite(R) & np.isfinite(T)
dau = ok & (R > 3 * rN) & (np.abs(P / v + 1) < 0.5)
print("daughter cells", dau.sum(), "T range", T[dau].min(), T[dau].max())
for T0 in np.arange(140, 262, 10):
    m = dau & (np.abs(T - T0) < 0.5)
    if m.sum() < 3:
        continue
    i = np.argmax(np.where(m, R, -1))
    a = R.flat[i]
    outer = m & (R > 0.5 * a)
    print("T=%5.1f cells %6d  a=max r=%10.4g  phi/v+1 at a: %+.2e  spread(phi/v)=%.1e  n r_N^2 at a: %+.3e spread(n r_N^2)=%.1e" % (
        T0, m.sum(), a, P.flat[i] / v + 1, np.std(P[outer] / v), n.flat[i] * rN**2, np.std(n[outer]) * rN**2))
