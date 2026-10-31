"""Round 11 (G11-01): strength of the Cauchy-horizon singularity of perturbed model-U holes.

Along an outgoing radial null geodesic u = const inside the hole (tangent k = dv/dlambda d_v, affine parameter dlambda = e^{2 sigma} dv
up to a constant) the only non-vanishing Ricci contraction is R_ab k^a k^b = R_vv (dv/dlambda)^2 with R_vv = -(D-2) nabla_v nabla_v r / r
= (D-2) phi_v^2 / h'(psi) (vv constraint).  Then
    lambda(v) = int e^{2 sigma} dv                                         (finite affine distance to the Cauchy horizon?)
    I_K(v)    = int R_kk dlambda = (D-2) int phi_v^2 e^{-2 sigma} / h' dv     (Krolak integral)
    I_T(v)    = int dlambda int R_kk dlambda'                              (Tipler integral; tidal distortion ~ I_T)
A singularity is Tipler-strong if I_T diverges as lambda -> lambda_CH, Krolak-strong if I_K diverges (Tipler 1977; Krolak 1987; for
null geodesics see Clarke & Krolak 1985).  A C^0 metric extension in a Kruskal-type coordinate V = -exp(-kappa v) exists if r stays
finite and 2 sigma + kappa v tends to a finite limit (Ori 1991 for Reissner-Nordstrom).
"""

from __future__ import annotations

import numpy as np

from . import doublenull10 as N


def ray_integrals(res: dict, D: int, frac=0.7, u_value=None):
    u, v = res["u"], res["v"]
    if u_value is None:
        ins = np.where(u > res["u_eh"] + 0.3)[0]
        i = ins[int(frac * (len(ins) - 1))]
    else:
        i = int(np.argmin(abs(u - u_value)))
    s, p, psi, r = res["s"][i], res["p"][i], res["psi"][i], res["r"][i]
    ok = np.isfinite(s) & np.isfinite(p) & np.isfinite(psi) & np.isfinite(r)
    v, s, p, psi, r = v[ok], s[ok], p[ok], psi[ok], r[ok]
    pv = np.gradient(p, v)
    h1 = 1.0 / np.clip(1.0 - psi, 1e-14, None) ** 2 if res["kind"] == "hayward" else (
        1.0 / np.clip(1.0 - psi, 1e-14, None) if res["kind"] == "dymnikova" else np.ones_like(psi))
    E2 = np.exp(2 * s)
    cum = lambda y: np.concatenate([[0.0], np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(v))])
    lam = cum(E2)
    IK = cum((D - 2) * pv**2 / h1 / E2)
    IT = cum(IK * E2)
    # Kruskal-type continuity: local slope of 2 sigma and the residual 2 sigma + kappa v with kappa from the late slope
    slope = np.gradient(2 * s, v)
    return dict(u=float(u[i]), v=v.tolist(), lam=lam.tolist(), I_K=IK.tolist(), I_T=IT.tolist(), sigma=s.tolist(), r=r.tolist(),
                slope_2sigma=slope.tolist())


def evolve(D, mu0, kind, amp, vmax, nu, nv, **kw):
    r1 = N.run_hyb(D, mu0, kind, nu=800, nv=3000, v_max=vmax, amp=amp, n_out_u=800, n_out_v=30, **kw)
    rend = r1["r"][:, -1]
    trapped = np.where(~(np.isfinite(rend) & (rend > r1["r0"])))[0]          # first row that does not escape
    k = int(trapped[0])
    ue = 0.5 * (r1["u"][k - 1] + r1["u"][k])
    res = N.run_hyb(D, mu0, kind, nu=nu, nv=nv, v_max=vmax, amp=amp, u_eh=ue, n_cluster=800, n_out_u=240, n_out_v=520, **kw)
    res["u_eh"] = ue
    return res


def summary(ri, at=(14, 18, 22, 24)):
    v = np.array(ri["v"])
    out = {}
    for vv in at:
        j = int(np.argmin(abs(v - vv)))
        out[vv] = dict(lam=ri["lam"][j], I_K=ri["I_K"][j], I_T=ri["I_T"][j], sigma=ri["sigma"][j], r=ri["r"][j], slope=ri["slope_2sigma"][j])
    return out


def classify(D, mu0, kind="hayward", amp=0.02, vmax=20.0, nu=2000, nv=8000, frac=0.7):
    """Run a pulse into a hole of mass mu0 and classify the interior outgoing ray: 'singular' if K grows by > 1e4 over its early value
    while the ray stays trapped; 'regular core' if the ray crosses an inner apparent horizon (f > 0) with bounded K."""
    from . import unified6l as U
    rm, rp = U.static_horizons(mu0, D, 1.0, kind)
    r0 = 1.6 * rp
    umax = 2 * (r0 - max(rm + 0.05, 1.03))
    r1 = N.run_hyb(D, mu0, kind, nu=800, nv=3000, v_max=vmax, amp=amp, u_max=umax, n_out_u=800, n_out_v=30)
    rend = r1["r"][:, -1]
    k = int(np.where(~(np.isfinite(rend) & (rend > r1["r0"])))[0][0])
    ue = 0.5 * (r1["u"][k - 1] + r1["u"][k])
    res = N.run_hyb(D, mu0, kind, nu=nu, nv=nv, v_max=vmax, amp=amp, u_max=umax, u_eh=ue, n_cluster=800, n_out_u=240, n_out_v=400)
    res["u_eh"] = ue
    u, v = res["u"], res["v"]
    ins = np.where(u > ue + 0.15)[0]
    i = ins[int(frac * (len(ins) - 1))]
    K = res["K"][i]; f = 1 - res["r"][i] ** 2 * res["psi"][i]
    ok = np.isfinite(K)
    K0 = np.nanmedian(K[(v > 1) & (v < 4)])
    Kmax = float(np.nanmax(K[ok])) if ok.any() else np.nan
    crossed = bool(np.any(f[ok][v[ok] > 2] > 0))
    v_cross = float(v[ok][v[ok] > 2][np.argmax(f[ok][v[ok] > 2] > 0)]) if crossed else None
    return dict(D=D, mu0=mu0, psi_minus=float(1 / rm**2), r_minus=float(rm), K_early=float(K0), K_max=Kmax, growth=float(Kmax / K0),
                crossed_inner_AH=crossed, v_cross=v_cross, v_last=float(v[ok][-1]),
                verdict="regular core (inner AH crossed)" if crossed and Kmax < 1e3 * K0 else ("singular growth" if Kmax > 1e4 * K0 else "undecided"))


def refined_run(D, mu0, kind="hayward", amp=0.02, vmax=22.0, nu=4000, nv=16000, passes=3, n_cluster=1200, umax=None, **kw):
    """Locate the event horizon by successive clustered runs (escaping vs non-escaping rows at v_max), then evolve.  Rays that hug the
    horizon fall over a u-width ~ exp(-kappa_+ v); an error delta in the horizon location limits the evolution to v ~ ln(1/delta)/kappa_+."""
    from . import unified6l as U
    rm, rp = U.static_horizons(mu0, D, 1.0, kind)
    r0 = 1.6 * rp
    umax = 2 * (r0 - max(rm + 0.05, 1.03)) if umax is None else umax
    ue = None
    for k in range(passes):
        res = N.run_hyb(D, mu0, kind, nu=nu if k == passes - 1 else 1200, nv=nv if k == passes - 1 else 4000, v_max=vmax, amp=amp,
                        u_max=umax, u_eh=ue, n_cluster=n_cluster if ue is not None else 0, n_out_u=10**6, n_out_v=50 if k < passes - 1 else 400, **kw)
        if k < passes - 1:
            rend = res["r"][:, -1]
            esc = np.isfinite(rend) & (rend > r0)
            j = int(np.where(~esc)[0][0])
            ue = 0.5 * (res["u"][j - 1] + res["u"][j])
    res["u_eh"] = ue
    return res
