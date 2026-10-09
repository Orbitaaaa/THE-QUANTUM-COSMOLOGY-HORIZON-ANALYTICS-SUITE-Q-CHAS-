"""Photon temperature transfer functions Delta_l(k) for the torus covariance, [TOP] Sec. 8 and Sec. 14.

CAMB 2.0.4 with the Planck 2018 best fit ([TOP] Sec. 14, boosted accuracy, every l computed
because ``lSampleBoost = 50``) supplies Delta_l(q) on its own wavenumber grid
(``get_cmb_transfer_data``, source 0 = temperature). We keep l = 2..40, interpolate in k with a
cubic spline through CAMB's samples, and tabulate it on a fine grid in k chi_* that covers
0 < k chi_* <= 200 with step 0.02 (CAMB's own grid beyond that).

CAMB artifact and repair: CAMB 2.0.4 returns Delta_l(q) = 0 exactly below k chi_* ~ 14.2-14.5
for 18 <= l <= 27 (the cut moves with the k sampling), although the Sachs-Wolfe term is still
sizable there (Delta_18 jumps from 0 to 11.5 % of its peak; a natural spline through the step
rings at about 1 % of the peak just above it). CAMB's own C_l carry the same truncation. Below
x = k chi_* < l the ratio Delta_l / j_l(x) is smooth (checked on the untruncated l = 16, 17), so for
every l whose first non-zero sample exceeds 1e-3 of its peak we continue
Delta_l(x) = j_l(x) (c0 + c1 x), fitted to the first four non-zero samples. This changes C_18 by
+0.08 %, C_19 by +0.015 % and the others by < 3e-5, and the torus KL divergences of
:mod:`qchas.cmb.covariance` by -1.4 % (cubic 1.5 D) to +0.1 % (-0.4 % to -0.7 % at 1.2-1.3 D); the
raw CAMB array is kept as :attr:`TransferFunctions.delta_q_camb` and ``load_transfer(repair=False)``
uses it. The paper ([TOP] Sec. 14: "the isotropic integral reproduces CAMB's C_l to 1e-15")
evidently used the raw arrays: raw transfer functions reproduce its printed BCC 1.2 D KL (0.50)
and cubic 1.3 D Pinsker bound (0.25), the repaired ones give 0.494 and 0.244. The covariance
module always uses the repaired functions.

Units: :meth:`TransferFunctions.delta` returns Delta_l(k) T_CMB in uK, so that

    C_l = 4 pi int dln k  P_R(k) Delta_l(k)^2        [uK^2]

with P_R the dimensionless curvature spectrum. The unit factor is found empirically as the
ratio of CAMB's unlensed scalar C_l to the dimensionless integral; it equals T_CMB^2 (uK^2)
to the quadrature accuracy and is stored in :attr:`TransferFunctions.unit_factor_measured`.
Results are cached in ``data/cache/transfer.npz`` (``$QCHAS_CACHE_DIR`` overrides the directory)
keyed by the parameters; an unwritable or corrupt cache only costs a recomputation.
"""

from __future__ import annotations

import contextlib
import hashlib
import json
import os
import zipfile
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.special import spherical_jn

from ..constants import PLANCK2018
from .camb_model import primordial_power

CACHE_DIR = Path(os.environ.get("QCHAS_CACHE_DIR", Path(__file__).resolve().parents[2] / "data" / "cache"))

LMIN_TRANSFER = 2
LMAX_TRANSFER = 40
#: CAMB settings of :func:`qchas.cmb.camb_model.camb_results` (lmax, AccuracyBoost, lSampleBoost)
CAMB_SETTINGS = {"lmax": 200, "accuracy_boost": 2.0, "l_sample_boost": 50.0}
#: fine tabulation grid in x = k chi_*: log-spaced below 1, step 0.02 up to 200
X_FINE_MAX = 200.0
X_FINE_STEP = 0.02
#: spline nodes are CAMB samples with k chi_* below this
X_SPLINE_MAX = 2500.0


def cache_key() -> str:
    """Hash of everything that determines the transfer functions."""
    import camb

    payload = {"params": PLANCK2018, "camb": camb.__version__, "settings": CAMB_SETTINGS,
               "lmin": LMIN_TRANSFER, "lmax": LMAX_TRANSFER, "x_fine": [X_FINE_MAX, X_FINE_STEP], "v": 3}
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()[:16]


@dataclass
class TransferFunctions:
    """CAMB temperature transfer functions in uK on CAMB's grid plus a spline in k."""

    ells: np.ndarray  # (nl,) = 2..40
    q: np.ndarray  # (nq,) CAMB wavenumbers, Mpc^-1
    delta_q: np.ndarray  # (nl, nq) Delta_l(q) T_CMB, uK
    chi_star: float  # Mpc
    t_cmb_uk: float
    cl_camb: np.ndarray  # CAMB unlensed scalar TT C_l, uK^2, l = 0..lmax_camb
    unit_factor_measured: float  # median over 2 <= l <= 30 of C_l^CAMB / (4 pi int P Delta_CAMB^2)
    delta_q_camb: np.ndarray  # (nl, nq) raw CAMB Delta_l(q) T_CMB before the low-k repair, uK
    repaired_ells: np.ndarray  # multipoles whose low-k truncation was repaired
    _spline: CubicSpline = field(init=False, repr=False)
    _fine: tuple[np.ndarray, np.ndarray] = field(init=False, repr=False)

    def __post_init__(self) -> None:
        sel = self.q * self.chi_star <= X_SPLINE_MAX
        self._spline = CubicSpline(self.q[sel], self.delta_q[:, sel].T, axis=0, bc_type="natural")
        x0 = self.q[0] * self.chi_star
        x = np.concatenate([np.geomspace(x0, 1.0, 200, endpoint=False), np.arange(1.0, X_FINE_MAX + 1e-9, X_FINE_STEP)])
        k = x / self.chi_star
        self._fine = (k, self._spline(k).T)

    @property
    def lmax(self) -> int:
        return int(self.ells[-1])

    @property
    def k_fine(self) -> np.ndarray:
        """Fine wavenumber grid (Mpc^-1) covering 0 < k chi_* <= 200."""
        return self._fine[0]

    @property
    def delta_fine(self) -> np.ndarray:
        """(nl, n_fine) Delta_l on :attr:`k_fine`, uK."""
        return self._fine[1]

    def delta_table(self, k, lmin: int = 2, lmax: int | None = None) -> np.ndarray:
        """Delta_l(k) for all l in lmin..lmax, shape k.shape + (nl,), uK.

        Below CAMB's first sample the small-k law Delta_l ~ k^l is used; above the spline
        range CAMB's samples are interpolated linearly.
        """
        lmax = self.lmax if lmax is None else lmax
        if lmin < self.ells[0] or lmax > self.lmax:
            raise ValueError(f"transfer functions cover l = {self.ells[0]}..{self.lmax}")
        k = np.asarray(k, dtype=float)
        flat = np.abs(k.ravel())
        il = slice(lmin - int(self.ells[0]), lmax - int(self.ells[0]) + 1)
        ls = self.ells[il]
        out = np.zeros((flat.size, ls.size))
        q0, qs = self.q[0], self._spline.x[-1]
        mid = (flat >= q0) & (flat <= qs)
        out[mid] = self._spline(flat[mid])[:, il]
        low = (flat > 0) & (flat < q0)
        if low.any():
            out[low] = self.delta_q[il, 0][None, :] * (flat[low, None] / q0) ** ls[None, :]
        high = flat > qs
        if high.any():
            for j, row in enumerate(self.delta_q[il]):
                out[high, j] = np.interp(flat[high], self.q, row, right=0.0)
        out[~np.isfinite(flat)] = np.nan
        return out.reshape(k.shape + (ls.size,))

    def delta(self, l, k) -> np.ndarray:
        """Delta_l(k) in uK, vectorized: ``l`` and ``k`` broadcast against each other."""
        l_arr, k_arr = np.broadcast_arrays(np.asarray(l, dtype=int), np.asarray(k, dtype=float))
        table = self.delta_table(k_arr.ravel(), int(self.ells[0]), self.lmax)
        idx = l_arr.ravel() - int(self.ells[0])
        if idx.size and (idx.min() < 0 or idx.max() >= self.ells.size):
            raise ValueError(f"l outside {self.ells[0]}..{self.lmax}")
        return table[np.arange(idx.size), idx].reshape(l_arr.shape)

    def _quadrature_grid(self, kmin: float, kmax: float) -> tuple[np.ndarray, np.ndarray]:
        """Wavenumbers and Delta table on [kmin, kmax]: fine grid, then CAMB samples."""
        kf, df = self._fine
        tail = self.q > kf[-1]
        k = np.concatenate([kf, self.q[tail]])
        d = np.concatenate([df, self.delta_q[:, tail]], axis=1)
        sel = (k > kmin) & (k < kmax)
        k, d = k[sel], d[:, sel]
        ends = [e for e in (kmin, kmax) if self.q[0] <= e <= self.q[-1]]
        if ends:
            ke = np.asarray(ends)
            de = self.delta_table(ke).T
            k = np.concatenate([k, ke])
            d = np.concatenate([d, de], axis=1)
            order = np.argsort(k)
            k, d = k[order], d[:, order]
        return k, d

    def cl_isotropic(self, lmin: int = 2, lmax: int = 30, kmin: float = 0.0, kmax: float = np.inf,
                     power=None) -> np.ndarray:
        """4 pi int_{kmin}^{kmax} dln k P_R(k) Delta_l(k)^2 for l = lmin..lmax, uK^2.

        ``power`` is a callable P_R(k); default the Planck 2018 power law.
        """
        if lmin < self.ells[0] or lmax > self.lmax or lmin > lmax:
            raise ValueError(f"transfer functions cover l = {self.ells[0]}..{self.lmax}")
        k, d = self._quadrature_grid(kmin, kmax)
        if k.size < 2:
            return np.zeros(lmax - lmin + 1)
        pr = (primordial_power if power is None else power)(k)
        il = slice(lmin - int(self.ells[0]), lmax - int(self.ells[0]) + 1)
        return 4.0 * np.pi * np.trapezoid(pr * d[il] ** 2, np.log(k), axis=1)

    def dl_isotropic(self, lmin: int = 2, lmax: int = 30, **kw) -> np.ndarray:
        ell = np.arange(lmin, lmax + 1)
        return ell * (ell + 1) * self.cl_isotropic(lmin, lmax, **kw) / (2.0 * np.pi)


def repair_low_k_truncation(x: np.ndarray, delta: np.ndarray, ells: np.ndarray, threshold: float = 1e-3,
                            n_fit: int = 4) -> tuple[np.ndarray, list[int]]:
    """Continue Delta_l = j_l(x)(c0 + c1 x) below CAMB's truncation where the step exceeds ``threshold``."""
    out = delta.copy()
    fixed = []
    for i, l in enumerate(ells):
        nz = int(np.flatnonzero(delta[i] != 0)[0])
        if nz == 0 or abs(delta[i, nz]) < threshold * np.max(np.abs(delta[i])):
            continue
        sl = slice(nz, nz + n_fit)
        c1, c0 = np.polyfit(x[sl], delta[i, sl] / spherical_jn(int(l), x[sl]), 1)
        out[i, :nz] = spherical_jn(int(l), x[:nz]) * (c0 + c1 * x[:nz])
        fixed.append(int(l))
    return out, fixed


def _compute() -> dict:
    """Run CAMB (shared cached results of :mod:`camb_model`) and extract l = 2..40."""
    from .camb_model import camb_results

    pars, res = camb_results(CAMB_SETTINGS["lmax"], CAMB_SETTINGS["accuracy_boost"])
    td = res.get_cmb_transfer_data(tp="scalar")
    L = np.asarray(td.L)
    rows = [int(np.flatnonzero(L == l)[0]) for l in range(LMIN_TRANSFER, LMAX_TRANSFER + 1)]
    q = np.asarray(td.q, dtype=float)
    delta = np.asarray(td.delta_p_l_k[0][rows], dtype=float)  # temperature source
    zstar = res.get_derived_params()["zstar"]
    chi = float(res.comoving_radial_distance(zstar))
    spec = res.get_cmb_power_spectra(CMB_unit="muK", raw_cl=True, lmax=CAMB_SETTINGS["lmax"])
    cl = np.asarray(spec["unlensed_scalar"][:, 0])
    t0 = float(pars.TCMB) * 1e6
    pr = primordial_power(q)
    raw = 4.0 * np.pi * np.trapezoid(pr * delta**2, np.log(q), axis=1)
    ratio = cl[LMIN_TRANSFER:31] / raw[: 31 - LMIN_TRANSFER]
    ells = np.arange(LMIN_TRANSFER, LMAX_TRANSFER + 1)
    fixed, ells_fixed = repair_low_k_truncation(q * chi, delta, ells)
    return dict(ells=ells, q=q, delta_q=fixed * t0, chi_star=chi, t_cmb_uk=t0, cl_camb=cl,
                unit_factor_measured=float(np.median(ratio)), delta_q_camb=delta * t0,
                repaired_ells=np.asarray(ells_fixed, dtype=int))


def _build(data: dict, repair: bool) -> TransferFunctions:
    if not repair:
        data = dict(data, delta_q=data["delta_q_camb"], repaired_ells=np.zeros(0, dtype=int))
    return TransferFunctions(**data)


@lru_cache(maxsize=2)
def load_transfer(refresh: bool = False, repair: bool = True) -> TransferFunctions:
    """Transfer functions, from ``data/cache/transfer.npz`` when its key matches.

    ``repair=False`` returns CAMB's raw arrays (with its low-k truncation for 18 <= l <= 27).
    """
    path = CACHE_DIR / "transfer.npz"
    key = cache_key()
    stored = None if refresh else load_npz(path)
    if stored is not None and str(stored.pop("key", "")) == key:
        data = {n: (v.item() if v.ndim == 0 else v) for n, v in stored.items()}
        return _build(data, repair)
    data = _compute()
    save_npz_atomic(path, key=np.array(key), **data)
    return _build(data, repair)


def load_npz(path) -> dict | None:
    """Arrays of a cache file, or None when it is missing or unreadable."""
    try:
        with np.load(path, allow_pickle=False) as z:
            return {n: z[n] for n in z.files}
    except (OSError, ValueError, KeyError, EOFError, zipfile.BadZipFile):
        return None


def save_npz_atomic(path, **arrays) -> bool:
    """Write a cache file atomically (unique temporary name, then rename); False if not writable."""
    path = Path(path)
    tmp = path.with_name(f".{path.stem}.{os.getpid()}.tmp.npz")
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        np.savez(tmp, **arrays)
        os.replace(tmp, path)
        return True
    except OSError:
        with contextlib.suppress(OSError):
            tmp.unlink(missing_ok=True)
        return False


def chi_star() -> float:
    """Comoving distance to last scattering (Mpc) from the same CAMB run (paper: 13 870)."""
    return load_transfer().chi_star


def delta(l, k) -> np.ndarray:
    """Delta_l(k) in uK for broadcastable ``l`` (2..40) and ``k`` (Mpc^-1)."""
    return load_transfer().delta(l, k)
