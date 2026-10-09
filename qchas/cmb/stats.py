"""Large-angle CMB statistics on the real sky: spectra, S_1/2, quadrupole p-values, low-l axes.

[TOP] Sec. 8 compares LambdaCDM with Planck 2018 on two classic anomalies, the low quadrupole and
S_1/2 = int_{-1}^{1/2} C(theta)^2 dcos(theta) with C(theta) = sum_l (2l+1)/(4 pi) C_l P_l(cos theta)
over 2 <= l <= 30; [TOP] Sec. 14 ("Statistics") evaluates it by 400-point Gauss-Legendre
quadrature on 10^5 full-sky realizations. For infinite LambdaCDM a full-sky estimate is
C_l chi^2_{2l+1} / (2l+1), independently for each l, which is how realizations are drawn here.

[TOP] Sec. 9.1 (i) measures the quadrupole and octopole axes, the directions n maximizing the
angular-momentum dispersion sum_m m^2 |a_lm(n)|^2 (de Oliveira-Costa et al. 2004, PRD 69,
063516). That sum equals <a|(n.J)^2|a> = n^T M n with M_ij = Re <a|{J_i, J_j}|a>/2, so the axis
is the top eigenvector of M; :func:`multipole_axis` also runs the HEALPix grid search with the
real Wigner rotations of :mod:`qchas.cmb.harmonics` and checks the two agree.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from functools import cache, lru_cache
from pathlib import Path

import numpy as np

from ..data.fetch import DEFAULT_RAW_DIR, available
from ..data.registry import DATASETS
from . import harmonics as H
from .maps import SkyData, sky_map_export  # noqa: F401  (sky_map_export re-exported for the exporter)

LMIN, LMAX = 2, 30
N_GAUSS_LEGENDRE = 400
N_REALIZATIONS = 100_000
SEED = 20260927

#: [TOP] Sec. 8 (text and table) and Sec. 14 ("Statistics"), Planck 2018 Commander spectrum
PAPER_STATS = {
    "S12_planck": 6694.0,
    "S12_lcdm_median": 34576.0,
    "p_S12_planck": 0.076,
    "p_C2_planck": 0.047,
    "D2_planck": 226.0,
    "D2_lcdm": 1016.0,
    "S12_copi_range": [3766.0, 8938.0],
    "source": "[TOP] Sec. 8 and its model table; Sec. 14 (Statistics)",
}
#: [TOP] Sec. 9.1 (i), Planck 2018 Commander map, Galactic (l, b) in degrees
PAPER_AXES = {
    "quadrupole": [242.8, 49.0],
    "octopole": [239.5, 62.9],
    "angle_deg": 14.0,
    "p_isotropic": 0.0297,
    "p_lcdm_sims": 0.0296,
    "source": "[TOP] Sec. 9.1 (i)",
}


@dataclass
class Spectrum:
    """A published TT spectrum D_l = l(l+1)C_l/2pi in uK^2."""

    name: str
    ell: np.ndarray
    D: np.ndarray
    err: np.ndarray
    err_minus: np.ndarray | None = None
    err_plus: np.ndarray | None = None
    extra: dict = field(default_factory=dict)
    meta: dict = field(default_factory=dict)

    def select(self, lmin: int = LMIN, lmax: int = LMAX) -> Spectrum:
        k = (self.ell >= lmin) & (self.ell <= lmax)
        if not np.array_equal(self.ell[k], np.arange(lmin, lmax + 1)):
            raise ValueError(f"{self.name}: multipoles {lmin}..{lmax} are not all present")

        def cut(a):
            return None if a is None else np.asarray(a)[k]

        return replace(self, ell=self.ell[k], D=self.D[k], err=self.err[k], err_minus=cut(self.err_minus),
                       err_plus=cut(self.err_plus), extra={key: cut(v) for key, v in self.extra.items()})

    def cl_band(self, lmin: int = LMIN, lmax: int = LMAX) -> np.ndarray:
        s = self.select(lmin, lmax)
        return cl_from_dl(s.ell, s.D)


def cl_from_dl(ell, D) -> np.ndarray:
    ell = np.asarray(ell, dtype=float)
    return 2 * np.pi * np.asarray(D, dtype=float) / (ell * (ell + 1))


def dl_from_cl(ell, cl) -> np.ndarray:
    ell = np.asarray(ell, dtype=float)
    return ell * (ell + 1) * np.asarray(cl, dtype=float) / (2 * np.pi)


# --- loaders -------------------------------------------------------------------------------

def load_wmap_tt(raw_dir: Path | None = None) -> Spectrum:
    """WMAP 9-yr combined TT spectrum: l, D_l, error, measurement part, cosmic-variance part (uK^2).

    For l <= 32 the file gives the maximum-likelihood value of the marginalized distribution.
    """
    raw_dir = Path(raw_dir or DEFAULT_RAW_DIR)
    path = available("wmap_tt", raw_dir)
    if path is None:
        raise FileNotFoundError(f"WMAP TT spectrum not found in {raw_dir}; run `qchas fetch`")
    t = np.loadtxt(path, comments="#")
    ds = DATASETS["wmap_tt"]
    return Spectrum(
        name="wmap", ell=t[:, 0].astype(int), D=t[:, 1], err=t[:, 2],
        extra={"err_measurement": t[:, 3], "err_cosmic_variance": t[:, 4]},
        meta={"file": ds.filename, "citation": ds.citation, "label": "WMAP 9-yr TT (published)",
              "estimator": "maximum likelihood of the marginalized distribution for l <= 32 (Hinshaw et al. 2013)"},
    )


def load_planck_tt(raw_dir: Path | None = None) -> Spectrum | None:
    """Planck 2018 TT spectrum (columns l, D_l, -dD_l, +dD_l in uK^2), or None if absent."""
    raw_dir = Path(raw_dir or DEFAULT_RAW_DIR)
    path = available("planck_tt", raw_dir)
    if path is None:
        return None
    t = np.loadtxt(path, comments="#")
    ds = DATASETS["planck_tt"]
    lo, hi = np.abs(t[:, 2]), np.abs(t[:, 3])
    return Spectrum(
        name="planck", ell=t[:, 0].astype(int), D=t[:, 1], err=0.5 * (lo + hi), err_minus=lo, err_plus=hi,
        meta={"file": ds.filename, "citation": ds.citation, "label": "Planck 2018 TT (published)",
              "estimator": "Commander likelihood for l < 30 (Planck 2018 results V)"},
    )


def lcdm_cl_band(lmin: int = LMIN, lmax: int = LMAX) -> np.ndarray:
    """CAMB Planck 2018 best-fit (unlensed) C_l for l = lmin..lmax, uK^2."""
    from .camb_model import lcdm_cl

    return np.asarray(lcdm_cl(max(lmax, 64)))[lmin: lmax + 1]


# --- correlation function and S_1/2 ------------------------------------------------------

def legendre_table(x, lmax: int) -> np.ndarray:
    """P_l(x) for l = 0..lmax, shape (lmax + 1, *x.shape)."""
    x = np.asarray(x, dtype=float)
    out = np.empty((lmax + 1,) + x.shape)
    out[0] = 1.0
    if lmax >= 1:
        out[1] = x
    for l in range(1, lmax):
        out[l + 1] = ((2 * l + 1) * x * out[l] - l * out[l - 1]) / (l + 1)
    return out


def _band_ell(cl_band: np.ndarray, lmin: int) -> np.ndarray:
    return np.arange(lmin, lmin + np.shape(cl_band)[-1])


def correlation_function(cl_band, theta_deg, lmin: int = LMIN) -> np.ndarray:
    """C(theta) = sum_l (2l+1)/(4pi) C_l P_l(cos theta) over the band; shape (..., n_theta)."""
    cl_band = np.asarray(cl_band, dtype=float)
    ell = _band_ell(cl_band, lmin)
    P = legendre_table(np.cos(np.radians(np.atleast_1d(theta_deg))), int(ell[-1]))[ell]
    return (cl_band * (2 * ell + 1) / (4 * np.pi)) @ P


@lru_cache(maxsize=8)
def s12_kernel(lmin: int = LMIN, lmax: int = LMAX, n_gl: int = N_GAUSS_LEGENDRE) -> np.ndarray:
    """Matrix I with S_1/2 = C^T I C: I_ll' = (2l+1)(2l'+1)/(4pi)^2 int_{-1}^{1/2} P_l P_l' dx (Gauss-Legendre)."""
    t, w = np.polynomial.legendre.leggauss(n_gl)
    x, w = 0.75 * t - 0.25, 0.75 * w
    ell = np.arange(lmin, lmax + 1)
    f = legendre_table(x, lmax)[ell] * ((2 * ell + 1) / (4 * np.pi))[:, None]
    k = (f * w) @ f.T
    k.setflags(write=False)
    return k


def s_half(cl_band, lmin: int = LMIN, n_gl: int = N_GAUSS_LEGENDRE) -> np.ndarray:
    """S_1/2 in uK^4 for C_l given over lmin..lmax on the last axis (vectorized over leading axes)."""
    cl_band = np.asarray(cl_band, dtype=float)
    k = s12_kernel(lmin, lmin + cl_band.shape[-1] - 1, n_gl)
    return np.einsum("...i,ij,...j->...", cl_band, k, cl_band)


def s_half_direct(cl_band, lmin: int = LMIN, n_gl: int = N_GAUSS_LEGENDRE) -> float:
    """S_1/2 by Gauss-Legendre quadrature of C(theta)^2 itself (cross-check of :func:`s_half`)."""
    t, w = np.polynomial.legendre.leggauss(n_gl)
    x, w = 0.75 * t - 0.25, 0.75 * w
    c = correlation_function(cl_band, np.degrees(np.arccos(x)), lmin)
    return float(np.sum(w * c**2))


# --- LambdaCDM realizations and p-values ----------------------------------------------------

def lcdm_realizations(model_cl_band, lmin: int = LMIN, n: int = N_REALIZATIONS, seed: int = SEED) -> np.ndarray:
    """Full-sky estimates C_l chi^2_{2l+1}/(2l+1), shape (n, n_l)."""
    model = np.asarray(model_cl_band, dtype=float)
    dof = 2 * _band_ell(model, lmin) + 1
    rng = np.random.default_rng(seed)
    return model * rng.chisquare(dof, size=(n, dof.size)) / dof


def anomaly_pvalues(obs_cl_band, model_cl_band, lmin: int = LMIN, sims: np.ndarray | None = None,
                    s_sims: np.ndarray | None = None) -> dict:
    """P(C_2 <= obs) and P(S_1/2 <= obs) under the model (simulated; C_2 also analytic)."""
    from scipy.stats import chi2

    obs = np.asarray(obs_cl_band, dtype=float)
    model = np.asarray(model_cl_band, dtype=float)
    if sims is None:
        sims = lcdm_realizations(model, lmin)
    if s_sims is None:
        s_sims = s_half(sims, lmin)
    s_obs = float(s_half(obs, lmin))
    out = {"S12": s_obs, "p_S12": float(np.mean(s_sims <= s_obs))}
    if lmin == 2:
        out.update({
            "D2": float(dl_from_cl(2, obs[0])),
            "p_C2": float(np.mean(sims[:, 0] <= obs[0])),
            "p_C2_analytic": float(chi2.cdf(5 * obs[0] / model[0], 5)),
        })
    return out


# --- multipole axes ------------------------------------------------------------------------

@cache
def _angular_momentum(l: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """J_x, J_y, J_z in the |l m> basis, m = -l..l (Condon-Shortley)."""
    m = np.arange(-l, l)
    jp = np.diag(np.sqrt((l - m) * (l + m + 1.0)), -1).astype(complex)
    jm = jp.conj().T
    return (jp + jm) / 2, (jp - jm) / 2j, np.diag(np.arange(-l, l + 1)).astype(complex)


def _complex_block(a_l: np.ndarray, l: int) -> np.ndarray:
    """Complex healpy-convention a_lm for m = -l..l from the real block of :mod:`harmonics`."""
    import healpy as hp

    alm = H.real_to_healpy(np.asarray(a_l, dtype=float), l, l, alm_lmax=l)
    pos = np.array([alm[hp.Alm.getidx(l, l, m)] for m in range(l + 1)])
    neg = (-1.0) ** np.arange(1, l + 1) * np.conj(pos[1:])
    return np.concatenate([neg[::-1], pos])


def axis_matrix(a_l: np.ndarray, l: int) -> np.ndarray:
    """Symmetric 3x3 M with sum_m m^2 |a_lm(n)|^2 = n^T M n (Galactic Cartesian axes)."""
    a = _complex_block(a_l, l)
    J = _angular_momentum(l)
    M = np.empty((3, 3))
    for i in range(3):
        for j in range(3):
            M[i, j] = np.real(a.conj() @ ((J[i] @ J[j] + J[j] @ J[i]) / 2) @ a)
    return M


def dispersion_by_rotation(a_l: np.ndarray, l: int, vecs: np.ndarray) -> np.ndarray:
    """sum_m m^2 (a^r_lm)^2 in the frame whose z axis is each of ``vecs`` (real Wigner rotations)."""
    v = np.atleast_2d(np.asarray(vecs, dtype=float))
    theta = np.arccos(np.clip(v[:, 2] / np.linalg.norm(v, axis=1), -1, 1))
    phi = np.arctan2(v[:, 1], v[:, 0])
    m2 = np.arange(-l, l + 1) ** 2.0
    out = np.empty(v.shape[0])
    for i, (th, ph) in enumerate(zip(theta, phi)):
        rot = H.rotation_from_euler(0.0, -th, -ph)  # maps n to +z
        out[i] = m2 @ (H.real_rotation_block(l, rot) @ a_l) ** 2
    return out


def _north(v: np.ndarray) -> np.ndarray:
    return -v if v[2] < 0 else v


def multipole_axis(a_l: np.ndarray, l: int, nside_search: int = 64) -> dict:
    """Axis of multipole ``l`` from its real block: HEALPix grid search, refined exactly.

    Returns Galactic (l, b) with b >= 0 for the exact axis (top eigenvector of :func:`axis_matrix`)
    and for the best grid direction, the normalized dispersion sum m^2|a|^2 / (l^2 sum |a|^2),
    and the separation between the two estimates.
    """
    import healpy as hp

    npix = hp.nside2npix(nside_search)
    vec = np.array(hp.pix2vec(nside_search, np.arange(npix))).T
    vec = vec[vec[:, 2] >= 0]
    disp = dispersion_by_rotation(a_l, l, vec)
    k = int(np.argmax(disp))
    grid = _north(vec[k])
    w, v = np.linalg.eigh(axis_matrix(a_l, l))
    exact = _north(v[:, -1])
    norm = l**2 * float(np.sum(np.asarray(a_l) ** 2))
    lg, bg = H.vec_to_gal(grid)
    le, be = H.vec_to_gal(exact)
    return {
        "ell": l, "l": float(le), "b": float(be), "vector": exact.tolist(),
        "l_grid": float(lg), "b_grid": float(bg), "nside_search": nside_search,
        "grid_offset_deg": float(np.degrees(np.arccos(min(1.0, abs(float(grid @ exact)))))),
        "dispersion": float(w[-1] / norm), "dispersion_grid": float(disp[k] / norm),
    }


def axis_angle_deg(v1, v2) -> float:
    """Angle between two headless axes, in [0, 90] degrees."""
    a, b = np.asarray(v1, float), np.asarray(v2, float)
    c = abs(a @ b) / (np.linalg.norm(a) * np.linalg.norm(b))
    return float(np.degrees(np.arccos(min(1.0, c))))


def quadrupole_octopole_alignment(alm_real: np.ndarray, lmin: int = LMIN, lmax: int = LMAX,
                                  nside_search: int = 64, source: str | None = None) -> dict:
    """Quadrupole and octopole axes of a stacked real a_lm vector and their alignment."""
    sl = H.block_slices(lmin, lmax)
    q = multipole_axis(alm_real[sl[2]], 2, nside_search)
    o = multipole_axis(alm_real[sl[3]], 3, nside_search)
    ang = axis_angle_deg(q["vector"], o["vector"])
    return {
        "source": source, "quadrupole": q, "octopole": o, "angle_deg": ang,
        "p_isotropic": float(1.0 - np.cos(np.radians(ang))), "computed": True,
        "method": "max_n sum_m m^2 |a_lm(n)|^2 (de Oliveira-Costa et al. 2004)",
        "paper": dict(PAPER_AXES, map="Planck-PR3-Commander"),
    }


# --- export-ready dictionaries -------------------------------------------------------------

def _spectrum_dict(s: Spectrum | None, lmin: int, lmax: int) -> dict | None:
    if s is None:
        return None
    b = s.select(lmin, lmax)
    out = {"ell": b.ell.tolist(), "D": b.D.tolist(), "err": b.err.tolist()}
    if b.err_minus is not None:
        out.update({"err_minus": b.err_minus.tolist(), "err_plus": b.err_plus.tolist()})
    out.update({"label": s.meta.get("label"), "estimator": s.meta.get("estimator"),
                "citation": s.meta.get("citation"), "file": s.meta.get("file")})
    return out


def observed_dict(wmap: Spectrum | None, planck: Spectrum | None, lmin: int = LMIN, lmax: int = LMAX) -> dict:
    """``cmb.json -> observed``: published spectra over lmin..lmax (Planck null when absent)."""
    return {"wmap": _spectrum_dict(wmap, lmin, lmax), "planck": _spectrum_dict(planck, lmin, lmax)}


def _hist(x: np.ndarray, bins: int = 60) -> dict:
    hi = float(np.quantile(x, 0.995))
    counts, edges = np.histogram(x, bins=bins, range=(0.0, hi))
    return {"edges": edges.tolist(), "counts": counts.tolist(), "overflow": int(np.sum(x > hi))}


def low_l_statistics(sky: SkyData | None, wmap: Spectrum | None, planck: Spectrum | None, *,
                     lmin: int = LMIN, lmax: int = LMAX, model_cl_band: np.ndarray | None = None,
                     n_realizations: int = N_REALIZATIONS, seed: int = SEED, nside_search: int = 64) -> dict:
    """``cmb.json -> stats``: D_2, S_1/2 and their LambdaCDM p-values for every available input.

    Inputs are the published spectra (``wmap``, ``planck``) and our own map extraction (``sky``).
    ``primary`` is Planck when present, else WMAP; ``D2_obs`` and ``S12_obs`` refer to it.
    """
    if sky is None and wmap is None and planck is None:
        raise ValueError("no observed spectrum or map given")
    model = lcdm_cl_band(lmin, lmax) if model_cl_band is None else np.asarray(model_cl_band, float)
    sims = lcdm_realizations(model, lmin, n_realizations, seed)
    s_sims = s_half(sims, lmin)
    inputs: dict[str, np.ndarray | None] = {
        "planck": None if planck is None else planck.cl_band(lmin, lmax),
        "wmap": None if wmap is None else wmap.cl_band(lmin, lmax),
        "map": None if sky is None else np.asarray(sky.cl),
    }
    res = {k: (None if c is None else anomaly_pvalues(c, model, lmin, sims, s_sims)) for k, c in inputs.items()}

    def pick(key: str) -> dict:
        return {k: (None if r is None else r[key]) for k, r in res.items()}

    primary = "planck" if res["planck"] else ("wmap" if res["wmap"] else "map")
    theta = np.arange(0.0, 180.01, 2.0)
    ell = np.arange(lmin, lmax + 1)
    pc = PAPER_STATS
    out = {
        "lmin": lmin, "lmax": lmax, "n_realizations": n_realizations, "seed": seed, "computed": True,
        "estimator": "full-sky C_l, LambdaCDM realizations C_l chi2_{2l+1}/(2l+1); S_1/2 by 400-point Gauss-Legendre",
        "primary": primary,
        "D2_obs": res[primary]["D2"] if lmin == 2 else None,
        "S12_obs": res[primary]["S12"],
        "D2_lcdm": float(dl_from_cl(2, model[0])) if lmin == 2 else None,
        "S12_lcdm": float(s_half(model, lmin)),
        "S12_lcdm_median": float(np.median(s_sims)),
        "S12_lcdm_quantiles": {q: float(np.quantile(s_sims, float(q))) for q in ("0.05", "0.16", "0.5", "0.84", "0.95")},
        "S12_lcdm_hist": _hist(s_sims),
        "D2": pick("D2"), "S12": pick("S12"),
        "p_C2": pick("p_C2"), "p_C2_analytic": pick("p_C2_analytic"), "p_S12": pick("p_S12"),
        "labels": {
            "planck": "Planck 2018 published TT (Commander, l < 30)",
            "wmap": "WMAP 9-yr published TT (maximum likelihood, l <= 32)",
            "map": None if sky is None else f"our full-sky extraction from the {sky.source_name} map",
        },
        "C_theta": {
            "theta_deg": theta.tolist(),
            "lcdm": correlation_function(model, theta, lmin).tolist(),
            **{k: (None if c is None else correlation_function(c, theta, lmin).tolist()) for k, c in inputs.items()},
        },
        "paper": dict(pc),
        "paper_check": {
            "note": "our LambdaCDM realizations evaluated at the Planck values printed in [TOP]",
            "p_S12_at_paper_S12": float(np.mean(s_sims <= pc["S12_planck"])),
            "p_C2_at_paper_D2": float(np.mean(sims[:, 0] <= cl_from_dl(2, pc["D2_planck"]))) if lmin == 2 else None,
            "S12_lcdm_median": float(np.median(s_sims)), "D2_lcdm": float(dl_from_cl(2, model[0])),
        },
        "lcdm_D": dl_from_cl(ell, model).tolist(),
    }
    if sky is not None and lmin <= 3 <= lmax:
        out["axes"] = quadrupole_octopole_alignment(sky.alm_real, lmin, lmax, nside_search, sky.source_name)
    return out


def extraction_check(sky: SkyData, wmap: Spectrum | None, planck: Spectrum | None) -> dict:
    """Our map D_l against the published spectrum of the same mission (and the [TOP] Sec. 9 Planck row)."""
    from .maps import PAPER_PIPELINE, SOURCE_PLANCK

    ref = planck if sky.source_name == SOURCE_PLANCK else wmap
    out = {"source": sky.source_name, "ell": sky.ell.tolist(), "map_D": sky.D_ell.tolist(),
           "published_D": None, "published": None, "computed": True,
           "paper": {k: PAPER_PIPELINE[k] for k in ("ell", "D_extracted", "D_published", "source")}}
    if ref is not None:
        out["published"] = ref.meta.get("label")
        out["published_D"] = ref.select(sky.lmin, sky.lmax).D.tolist()
    return out


def cmb_fields(sky: SkyData | None, wmap: Spectrum | None = None, planck: Spectrum | None = None, *,
               lmin: int = LMIN, lmax: int = LMAX, **kw) -> dict:
    """The ``source``, ``ell``, ``observed``, ``map_D``, ``lcdm_D`` and ``stats`` fields of ``cmb.json``."""
    stats = low_l_statistics(sky, wmap, planck, lmin=lmin, lmax=lmax, **kw)
    if sky is not None:
        stats["extraction_check"] = extraction_check(sky, wmap, planck)
    return {
        "source": None if sky is None else sky.source_name,
        "ell": list(range(lmin, lmax + 1)),
        "observed": observed_dict(wmap, planck, lmin, lmax),
        "map_D": None if sky is None else np.asarray(sky.D_ell).tolist(),
        "lcdm_D": stats.pop("lcdm_D"),
        "stats": stats,
    }
