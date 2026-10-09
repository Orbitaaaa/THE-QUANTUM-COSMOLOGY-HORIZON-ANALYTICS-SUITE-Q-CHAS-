"""Real-sky CMB maps to low-l real multipoles, following the [TOP] Sec. 9 / Sec. 14 pipeline.

[TOP] Sec. 9 ("Data") and Sec. 14 ("Planck pipeline"): the inpainted Commander temperature is
averaged to Nside 64 in the NESTED scheme and the real multipoles 2 <= l <= 30 are extracted by
HEALPix quadrature with four Jacobi iterations up to l = 64 and a pixel-window correction
"measured by simulation"; on a simulated sky processed identically the residual is 0.26 % of
the signal.

Steps implemented here, for either the Planck PR3 Commander map (preferred, 5' FWHM) or the
WMAP 9-year ILC map (1 deg FWHM, Bennett et al. 2013, ApJS 208, 20):

1. read the full-resolution map in its native NESTED ordering, convert to uK;
2. average every block of (N_in/64)^2 NESTED children: the input pixels are pixel averages, so
   this is exactly the Nside-64 pixel average and the only pixel window left is the Nside-64 one;
3. remove the full-sky monopole and dipole;
4. ``healpy.map2alm`` with ``iter=4`` up to ``lmax_alm`` (64 by default, as in the paper);
5. divide by b_l p_l (input Gaussian beam times Nside-64 pixel window);
6. convert to the stacked real vector of :mod:`qchas.cmb.harmonics` (``healpy_to_real``).

Two choices differ in form from the paper's wording: p_l is the analytic HEALPix window (healpy's
file, or :func:`exact_pixel_window`), and :func:`simulate_pipeline` checks it instead of
measuring it (per-l gains within 3e-3 of 1); the beam division is ours, negligible for the 5'
Commander beam (b_30 = 0.9998) and essential for the 1 deg ILC (b_30 = 0.85).

The same steps run on Gaussian LambdaCDM skies (:func:`simulate_pipeline`) to measure the residual.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path

import numpy as np

from ..constants import PLANCK2018
from ..data.fetch import DEFAULT_RAW_DIR, available, load_provenance
from ..data.registry import DATASETS
from . import harmonics as H

#: WMAP 9-yr ILC map resolution: smoothed to a 1 deg FWHM Gaussian beam (Bennett et al. 2013, ApJS 208, 20)
WMAP_ILC_BEAM_FWHM_ARCMIN = 60.0
#: Planck PR3 Commander CMB map resolution: 5' FWHM Gaussian (Planck 2018 results IV, A&A 641, A4)
PLANCK_COMMANDER_BEAM_FWHM_ARCMIN = 5.0

SOURCE_PLANCK = "Planck-PR3-Commander"
SOURCE_WMAP = "WMAP9-ILC"

#: Pipeline numbers printed in [TOP] Sec. 9 for the Planck 2018 Commander map
PAPER_PIPELINE = {
    "residual_fraction": 0.0026,
    "ell": [2, 3, 4, 5, 6],
    "D_extracted": [213.0, 898.0, 704.0, 1411.0, 570.0],
    "D_published": [226.0, 937.0, 692.0, 1502.0, 558.0],
    "source": "[TOP] Sec. 9",
}

_UNSEEN_LIMIT = -1e30


@dataclass
class FullResMap:
    """A full-resolution temperature map in NESTED ordering, in uK."""

    name: str
    map_uK: np.ndarray
    nside: int
    beam_fwhm_arcmin: float
    meta: dict = field(default_factory=dict)


@dataclass
class SkyData:
    """Low-l sky at ``nside`` (RING, Galactic, uK) with its real multipoles lmin..lmax."""

    source_name: str
    nside: int
    lmin: int
    lmax: int
    map_uK: np.ndarray
    alm_real: np.ndarray
    ell: np.ndarray
    cl: np.ndarray
    D_ell: np.ndarray
    lowl_map_uK: np.ndarray
    mask: np.ndarray
    mask_fraction: np.ndarray | None
    fsky: float
    meta: dict = field(default_factory=dict)


@dataclass
class Extraction:
    """Output of :func:`extract_low_l`."""

    map_uK: np.ndarray
    alm: np.ndarray
    alm_lmax: int
    alm_real: np.ndarray
    ell: np.ndarray
    cl: np.ndarray
    monopole_uK: float
    dipole_uK: np.ndarray
    info: dict


# --- small helpers -----------------------------------------------------------------------

def _unit_factor(unit: str | None) -> float:
    """Multiplier from a thermodynamic (CMB) FITS temperature unit to uK."""
    u = (unit or "").strip().lower().replace("_", " ").replace(",", " ")
    if any(w in u.split() for w in ("rj", "antenna", "brightness")) or "mjy" in u:
        raise ValueError(f"temperature unit {unit!r} is not thermodynamic (CMB) temperature")
    head = u.split()[0] if u.split() else ""
    if head in ("uk", "muk", "microk"):
        return 1.0
    if head == "mk":
        return 1e3
    if head in ("k", "kcmb") or head.startswith("k"):
        return 1e6
    raise ValueError(f"unrecognised temperature unit {unit!r}")


def _header(path: Path, hdu: int = 1) -> dict:
    from astropy.io import fits

    h = fits.getheader(path, hdu)
    cols = [h.get(f"TTYPE{i}", "").strip() for i in range(1, int(h.get("TFIELDS", 0)) + 1)]
    units = [h.get(f"TUNIT{i}", "").strip() for i in range(1, len(cols) + 1)]
    return {
        "ordering": str(h.get("ORDERING", "RING")).strip().upper(),
        "nside": int(h["NSIDE"]),
        "coordsys": str(h.get("COORDSYS", "")).strip() or None,
        "columns": cols,
        "units": units,
    }


def _provenance(key: str, raw_dir: Path) -> dict | None:
    rec = load_provenance(raw_dir).get(key)
    if rec is None:
        return None
    return {k: rec.get(k) for k in ("filename", "url", "bytes", "sha256", "retrieved_utc")}


def gaussian_beam(fwhm_arcmin: float, lmax: int) -> np.ndarray:
    """Gaussian beam transfer b_l = exp(-l(l+1) sigma^2 / 2), l = 0..lmax."""
    sigma = np.radians(fwhm_arcmin / 60.0) / np.sqrt(8.0 * np.log(2.0))
    ell = np.arange(lmax + 1)
    return np.exp(-0.5 * ell * (ell + 1) * sigma**2)


# --- pixel window --------------------------------------------------------------------------

def _legendre_mean(x: np.ndarray, lmax: int) -> np.ndarray:
    """Row means of P_l(x) for l = 0..lmax; ``x`` has shape (n_pix, n_sub)."""
    out = np.empty((lmax + 1, x.shape[0]))
    p0, p1 = np.ones_like(x), x.copy()
    out[0] = 1.0
    if lmax >= 1:
        out[1] = p1.mean(axis=1)
    for l in range(1, lmax):
        p0, p1 = p1, ((2 * l + 1) * x * p1 - l * p0) / (l + 1)
        out[l + 1] = p1.mean(axis=1)
    return out


def _subpixel_window(nside: int, lmax: int, sub: int) -> np.ndarray:
    import healpy as hp

    npix = hp.nside2npix(nside)
    theta, phi = hp.pix2ang(nside, np.arange(npix))
    # one octant (north, 0 <= phi < pi/2) represents all pixels by the HEALPix symmetries;
    # pixels centred on the equator are shared with the south and count one half
    north = np.where((theta < np.pi / 2 - 1e-12) & (phi < np.pi / 2))[0]
    equator = np.where((np.abs(theta - np.pi / 2) <= 1e-12) & (phi < np.pi / 2))[0]
    pix = hp.ring2nest(nside, np.concatenate([north, equator]))
    weight = np.concatenate([np.ones(north.size), np.full(equator.size, 0.5)])
    k = sub * sub
    acc = np.zeros(lmax + 1)
    for chunk in np.array_split(np.arange(pix.size), max(1, pix.size * k // 400_000)):
        p = pix[chunk]
        centre = np.array(hp.pix2vec(nside, p, nest=True)).T
        child = (p[:, None] * k + np.arange(k)[None, :]).ravel()
        v = np.array(hp.pix2vec(nside * sub, child, nest=True)).T.reshape(p.size, k, 3)
        x = np.clip(np.einsum("pkj,pj->pk", v, centre), -1.0, 1.0)
        acc += _legendre_mean(x, lmax) @ weight[chunk]
    return acc / weight.sum()


@lru_cache(maxsize=8)
def exact_pixel_window(nside: int, lmax: int) -> np.ndarray:
    """Linear-response pixel window W_l = < mean_{x in p} P_l(n_p . x) >_p, l = 0..lmax.

    This is the factor by which pixel averaging multiplies a_lm on average (the quantity a
    simulation measures). The pixel integral is done on the 16x16 and 8x8 NESTED children and
    Richardson-extrapolated; it agrees with the healpy pixel-window file for Nside 64 to 1e-5 at
    l <= 30.
    """
    w16 = _subpixel_window(nside, lmax, 16)
    w8 = _subpixel_window(nside, lmax, 8)
    w = (4.0 * w16 - w8) / 3.0
    w.setflags(write=False)
    return w


@lru_cache(maxsize=16)
def _healpy_pixwin(nside: int, lmax: int) -> np.ndarray | None:
    """``healpy.pixwin`` or None when its data file is neither cached nor downloadable."""
    import healpy as hp

    try:
        w = np.asarray(hp.pixwin(nside, lmax=lmax), dtype=float)
    except (OSError, ValueError):  # urllib errors are OSError subclasses
        return None
    w.setflags(write=False)
    return w


def _pixel_window_cached(nside: int, lmax: int) -> tuple[np.ndarray, str]:
    w = _healpy_pixwin(nside, lmax)
    return (w, "healpy") if w is not None else (exact_pixel_window(nside, lmax), "exact")


def pixel_window(nside: int, lmax: int) -> tuple[np.ndarray, str]:
    """Temperature pixel window p_l (l = 0..lmax) and the method used ('healpy' or 'exact').

    ``healpy.pixwin`` is used when its data file is available; otherwise the window is computed
    by :func:`exact_pixel_window`.
    """
    w, method = _pixel_window_cached(int(nside), int(lmax))
    return w.copy(), method


# --- loaders -------------------------------------------------------------------------------

def _read_nested(path: Path, field_index: int, memmap: bool = False) -> np.ndarray:
    import healpy as hp

    return hp.read_map(str(path), field=field_index, nest=True, memmap=memmap)


def load_wmap_ilc(raw_dir: Path | None = None) -> FullResMap:
    """WMAP 9-yr ILC map (Nside 512) in NESTED ordering and uK, beam 1 deg FWHM.

    The file is NESTED with unit 'mK, thermodynamic'; both are read from its header.
    """
    raw_dir = Path(raw_dir or DEFAULT_RAW_DIR)
    path = available("wmap_ilc", raw_dir)
    if path is None:
        raise FileNotFoundError(f"WMAP ILC map not found in {raw_dir}; run `qchas fetch`")
    hdr = _header(path)
    col = hdr["columns"].index("TEMPERATURE") if "TEMPERATURE" in hdr["columns"] else 0
    factor = _unit_factor((hdr["units"][col] if hdr["units"] else "") or "mK")
    m = np.asarray(_read_nested(path, col), dtype=np.float64) * factor
    ds = DATASETS["wmap_ilc"]
    return FullResMap(
        name=SOURCE_WMAP, map_uK=m, nside=hdr["nside"], beam_fwhm_arcmin=WMAP_ILC_BEAM_FWHM_ARCMIN,
        meta={
            "file": ds.filename, "citation": ds.citation, "column": hdr["columns"][col],
            "ordering_in": hdr["ordering"], "coordsys": hdr["coordsys"] or "G (WMAP Galactic)",
            "unit_in": hdr["units"][col] if hdr["units"] else None, "unit_factor_to_uK": factor,
            "beam_fwhm_arcmin": WMAP_ILC_BEAM_FWHM_ARCMIN,
            "beam_source": "Bennett et al. 2013, ApJS 208, 20: 9-yr ILC map smoothed to 1 deg FWHM",
            "provenance": _provenance("wmap_ilc", raw_dir),
        },
    )


def load_planck_commander(raw_dir: Path | None = None, memmap: bool = True) -> FullResMap | None:
    """Planck PR3 Commander temperature (Nside 2048, NESTED, uK, 5' FWHM), or None if absent.

    The inpainted column ``I_STOKES_INP`` is used when the file has it, else field 0
    (``I_STOKES``); the choice is recorded in ``meta['column']``.
    """
    raw_dir = Path(raw_dir or DEFAULT_RAW_DIR)
    path = available("planck_commander", raw_dir)
    if path is None:
        return None
    hdr = _header(path)
    cols = [c.upper() for c in hdr["columns"]]
    col = cols.index("I_STOKES_INP") if "I_STOKES_INP" in cols else 0
    unit = (hdr["units"][col] if hdr["units"] else "") or "K_CMB"
    factor = _unit_factor(unit)
    m = np.asarray(_read_nested(path, col, memmap=memmap), dtype=np.float32) * np.float32(factor)
    ds = DATASETS["planck_commander"]
    return FullResMap(
        name=SOURCE_PLANCK, map_uK=m, nside=hdr["nside"], beam_fwhm_arcmin=PLANCK_COMMANDER_BEAM_FWHM_ARCMIN,
        meta={
            "file": ds.filename, "citation": ds.citation, "column": hdr["columns"][col],
            "ordering_in": hdr["ordering"], "coordsys": hdr["coordsys"],
            "unit_in": unit, "unit_factor_to_uK": factor,
            "beam_fwhm_arcmin": PLANCK_COMMANDER_BEAM_FWHM_ARCMIN,
            "beam_source": "Planck 2018 results IV (A&A 641, A4): CMB maps at 5 arcmin FWHM",
            "provenance": _provenance("planck_commander", raw_dir),
        },
    )


def load_kq75_mask(nside: int = 64, threshold: float = 0.5,
                   raw_dir: Path | None = None) -> tuple[np.ndarray, np.ndarray, float] | None:
    """WMAP KQ75 mask at ``nside`` (RING): (uint8 mask, kept fraction per pixel, native f_sky).

    A pixel is kept (1) when at least ``threshold`` of its native pixels are kept. None if absent.
    """
    raw_dir = Path(raw_dir or DEFAULT_RAW_DIR)
    path = available("wmap_kq75", raw_dir)
    if path is None:
        return None
    m = np.asarray(_read_nested(path, 0), dtype=np.float64)
    native = (m > 0.5).astype(np.float64)
    frac = degrade_nested(native, nside)
    return (frac >= threshold).astype(np.uint8), frac, float(native.mean())


# --- pipeline ------------------------------------------------------------------------------

def degrade_nested(map_nested: np.ndarray, nside_out: int) -> np.ndarray:
    """Average a NESTED map to ``nside_out`` and return it in RING ordering.

    UNSEEN or non-finite children are ignored; an output pixel with no valid child is UNSEEN.
    """
    import healpy as hp

    npix_out = 12 * nside_out**2
    m = np.asarray(map_nested)
    if m.size % npix_out:
        raise ValueError("input map size is not a NESTED refinement of nside_out")
    blocks = m.reshape(npix_out, m.size // npix_out)
    bad = ~np.isfinite(blocks) | (blocks < _UNSEEN_LIMIT)
    if not bad.any():
        out = blocks.mean(axis=1, dtype=np.float64)
    else:
        good = ~bad
        n = good.sum(axis=1)
        s = np.where(good, blocks, 0.0).sum(axis=1, dtype=np.float64)
        out = np.where(n > 0, s / np.maximum(n, 1), hp.UNSEEN)
    return hp.reorder(out, n2r=True)


def deconvolved_alm(map_ring: np.ndarray, beam_fwhm_arcmin: float, lmax_alm: int,
                    iterations: int = 4, input_window: np.ndarray | None = None) -> tuple[np.ndarray, str]:
    """``map2alm`` (Jacobi ``iterations``) divided by beam x pixel window; returns (alm, pixwin method).

    ``input_window`` is an optional extra transfer function (l = 0..lmax_alm) to divide out.
    """
    import healpy as hp

    nside = hp.npix2nside(map_ring.size)
    alm = hp.map2alm(np.asarray(map_ring, dtype=np.float64), lmax=lmax_alm, iter=iterations)
    pw, method = pixel_window(nside, lmax_alm)
    transfer = gaussian_beam(beam_fwhm_arcmin, lmax_alm) * pw[: lmax_alm + 1]
    if input_window is not None:
        transfer = transfer * np.asarray(input_window)[: lmax_alm + 1]
    return hp.almxfl(alm, 1.0 / transfer), method


def cl_from_real(alm_real: np.ndarray, lmin: int, lmax: int) -> np.ndarray:
    """Full-sky C_l = sum_m (a^r_lm)^2 / (2l+1) for l = lmin..lmax."""
    return np.array([np.mean(alm_real[s] ** 2) for s in H.block_slices(lmin, lmax).values()])


def dl_from_cl(ell: np.ndarray, cl: np.ndarray) -> np.ndarray:
    ell = np.asarray(ell, dtype=float)
    return ell * (ell + 1) * np.asarray(cl) / (2 * np.pi)


def extract_low_l(map_nested: np.ndarray, beam_fwhm_arcmin: float, *, nside: int = 64, lmin: int = 2,
                  lmax: int = 30, lmax_alm: int | None = None, iterations: int = 4,
                  input_window: np.ndarray | None = None) -> Extraction:
    """Run steps 2-6 of the module pipeline on a full-resolution NESTED map in uK."""
    import healpy as hp

    lmax_alm = max(lmax, nside) if lmax_alm is None else int(lmax_alm)
    if lmax_alm < lmax:
        raise ValueError("lmax_alm must be >= lmax")
    m = degrade_nested(map_nested, nside)
    n_unseen = int(np.sum(m < _UNSEEN_LIMIT))
    if n_unseen:
        m = np.where(m < _UNSEEN_LIMIT, 0.0, m)
    m, mono, dip = hp.remove_dipole(m, fitval=True)
    alm, method = deconvolved_alm(m, beam_fwhm_arcmin, lmax_alm, iterations, input_window)
    vec = H.healpy_to_real(alm, lmin, lmax, lmax_alm)
    ell = np.arange(lmin, lmax + 1)
    info = {
        "nside": nside, "lmax_alm": lmax_alm, "iterations": iterations,
        "degrade": f"NESTED block average to Nside {nside} (exact pixel average of a pixel-averaged input)",
        "pixel_window": {"nside": nside, "method": method},
        "beam_fwhm_arcmin": beam_fwhm_arcmin, "unseen_pixels_zeroed": n_unseen,
        "extra_window_divided": input_window is not None,
    }
    return Extraction(map_uK=np.asarray(m, dtype=np.float64), alm=alm, alm_lmax=lmax_alm, alm_real=vec,
                      ell=ell, cl=cl_from_real(vec, lmin, lmax), monopole_uK=float(mono),
                      dipole_uK=np.asarray(dip, dtype=float), info=info)


def lowl_map(alm_real: np.ndarray, lmin: int, lmax: int, nside: int) -> np.ndarray:
    """Synthesis of the real multipoles lmin..lmax at ``nside`` (RING, no beam or pixel window)."""
    import healpy as hp

    return hp.alm2map(H.real_to_healpy(alm_real, lmin, lmax), nside, lmax=lmax)


def prepare_sky(source: str = "auto", nside: int = 64, lmin: int = 2, lmax: int = 30, *,
                raw_dir: Path | None = None, lmax_alm: int | None = None, iterations: int = 4,
                mask_threshold: float = 0.5) -> SkyData:
    """Low-l sky from Planck PR3 Commander (preferred when on disk) or the WMAP 9-yr ILC.

    ``source`` is 'auto', 'planck' or 'wmap'. The returned map is Nside ``nside`` RING, Galactic,
    uK with monopole and dipole removed; ``alm_real`` uses the :mod:`qchas.cmb.harmonics` convention.
    Without the KQ75 file ``mask`` is all ones, ``mask_fraction`` is None and a warning is recorded.
    """
    raw_dir = Path(raw_dir or DEFAULT_RAW_DIR)
    src = source.lower()
    warnings: list[str] = []
    full: FullResMap | None = None
    if src in ("auto", "planck"):
        full = load_planck_commander(raw_dir)
        if full is None:
            if src == "planck":
                raise FileNotFoundError(f"Planck Commander map not found in {raw_dir}")
            warnings.append("Planck PR3 Commander map not available; using the WMAP 9-year ILC map")
    if full is None:
        if src not in ("auto", "wmap"):
            raise ValueError(f"unknown source {source!r}")
        full = load_wmap_ilc(raw_dir)

    ex = extract_low_l(full.map_uK, full.beam_fwhm_arcmin, nside=nside, lmin=lmin, lmax=lmax,
                       lmax_alm=lmax_alm, iterations=iterations)
    name, nside_in, meta = full.name, full.nside, dict(full.meta)
    del full
    if name == SOURCE_PLANCK and str(meta.get("column", "")).upper() != "I_STOKES_INP":
        warnings.append(f"Planck Commander file has no I_STOKES_INP column; using {meta.get('column')} "
                        "(not the inpainted map of [TOP] Sec. 9)")
    if ex.info["unseen_pixels_zeroed"]:
        warnings.append(f"{name}: {ex.info['unseen_pixels_zeroed']} Nside-{nside} pixels had no valid input and "
                        "were set to zero before the full-sky harmonic analysis")

    mask_out = load_kq75_mask(nside, mask_threshold, raw_dir)
    if mask_out is None:
        mask, frac, fsky_native = np.ones(12 * nside**2, dtype=np.uint8), None, None
        warnings.append("WMAP KQ75 mask not available; mask set to the full sky")
    else:
        mask, frac, fsky_native = mask_out
    meta.update({
        "source": name, "nside_in": nside_in, "nside": nside, "ordering": "RING", "frame": "galactic",
        "unit": "uK", "lmin": lmin, "lmax": lmax, "extraction": ex.info,
        "monopole_removed_uK": ex.monopole_uK, "dipole_removed_uK": ex.dipole_uK.tolist(),
        "mask": {"name": "WMAP KQ75 (9-yr)", "threshold": mask_threshold, "fsky_native": fsky_native,
                 "fsky": float(mask.mean())} if frac is not None else None,
        "warnings": warnings,
    })
    return SkyData(
        source_name=name, nside=nside, lmin=lmin, lmax=lmax, map_uK=ex.map_uK, alm_real=ex.alm_real,
        ell=ex.ell, cl=ex.cl, D_ell=dl_from_cl(ex.ell, ex.cl), lowl_map_uK=lowl_map(ex.alm_real, lmin, lmax, nside),
        mask=mask, mask_fraction=frac, fsky=float(mask.mean()), meta=meta,
    )


# --- validation by simulation --------------------------------------------------------------

@lru_cache(maxsize=4)
def _camb_cl_standard(lmax: int) -> np.ndarray:
    import camb

    p = PLANCK2018
    pars = camb.set_params(H0=p["H0"], ombh2=p["ombh2"], omch2=p["omch2"], tau=p["tau"], As=p["As"],
                           ns=p["ns"], mnu=p["mnu"], lmax=lmax)
    spec = camb.get_results(pars).get_cmb_power_spectra(CMB_unit="muK", raw_cl=True, lmax=lmax)
    return np.asarray(spec["unlensed_scalar"][: lmax + 1, 0])


def lcdm_cl_simulation(lmax: int, boosted_lmax: int = 200) -> np.ndarray:
    """Unlensed LambdaCDM TT C_l (uK^2, l = 0..lmax) for simulated skies.

    l <= ``boosted_lmax`` comes from :func:`qchas.cmb.camb_model.lcdm_cl` (paper accuracy); the
    tail from a standard-accuracy CAMB run with the same Planck 2018 parameters (the boosted
    settings need tens of GB beyond l ~ 2000). The two agree to 0.1 % for l <= 1500.
    """
    cl = _camb_cl_standard(max(int(lmax), 2)).copy()[: lmax + 1]
    if boosted_lmax >= 2:
        from .camb_model import lcdm_cl

        lb = min(int(boosted_lmax), int(lmax))
        cl[: lb + 1] = lcdm_cl(lb)
    return cl


def gaussian_alm(cl: np.ndarray, lmax: int, rng: np.random.Generator, lmin: int = 2) -> np.ndarray:
    """Isotropic Gaussian healpy alm (l <= lmax) with power ``cl``; multipoles below lmin are zero."""
    import healpy as hp

    ell, m = hp.Alm.getlm(lmax)
    s = np.sqrt(np.asarray(cl, dtype=float)[ell])
    re, im = rng.standard_normal(ell.size), rng.standard_normal(ell.size)
    alm = np.where(m == 0, s * re, s * (re + 1j * im) / np.sqrt(2.0))
    alm[ell < lmin] = 0.0
    return alm


def simulate_pipeline(nside_in: int = 512, beam_fwhm_arcmin: float = WMAP_ILC_BEAM_FWHM_ARCMIN, *,
                      seed: int = 0, nside: int = 64, lmin: int = 2, lmax: int = 30,
                      lmax_alm: int | None = None, iterations: int = 4, cl: np.ndarray | None = None,
                      boosted_lmax: int = 200) -> dict:
    """Process a Gaussian LambdaCDM sky through the identical pipeline and compare a_lm.

    The sky is drawn to l = 3 nside_in - 1, convolved with the beam and, when the healpy
    pixel-window file for ``nside_in`` is available, with that pixel window (so that its pixels
    approximate pixel averages like a real map); otherwise it is point-sampled. Returns the
    relative residual ||a_out - a_in|| / ||a_in|| over lmin..lmax, per-l residuals and gains.
    """
    import healpy as hp

    lmax_in = 3 * nside_in - 1
    if cl is None:
        cl = lcdm_cl_simulation(lmax_in, boosted_lmax)
    rng = np.random.default_rng(seed)
    alm_in = gaussian_alm(cl, lmax_in, rng, lmin=lmin)
    transfer = gaussian_beam(beam_fwhm_arcmin, lmax_in)
    pw_in = _healpy_pixwin(nside_in, lmax_in)
    if pw_in is not None:
        transfer = transfer * pw_in[: lmax_in + 1]
    sky = hp.alm2map(hp.almxfl(alm_in, transfer), nside_in, lmax=lmax_in)
    sky = hp.reorder(sky, r2n=True)
    ex = extract_low_l(sky, beam_fwhm_arcmin, nside=nside, lmin=lmin, lmax=lmax, lmax_alm=lmax_alm,
                       iterations=iterations)
    a_in = H.healpy_to_real(alm_in, lmin, lmax, lmax_in)
    a_out = ex.alm_real
    per_l, gain = [], []
    for s in H.block_slices(lmin, lmax).values():
        per_l.append(float(np.linalg.norm(a_out[s] - a_in[s]) / np.linalg.norm(a_in[s])))
        gain.append(float(a_out[s] @ a_in[s] / (a_in[s] @ a_in[s])))
    return {
        "nside_in": nside_in, "beam_fwhm_arcmin": beam_fwhm_arcmin, "seed": seed, "nside": nside,
        "lmin": lmin, "lmax": lmax, "lmax_alm": ex.alm_lmax, "iterations": iterations,
        "input_pixwin": pw_in is not None, "pixel_window_method": ex.info["pixel_window"]["method"],
        "relative_residual": float(np.linalg.norm(a_out - a_in) / np.linalg.norm(a_in)),
        "ell": list(range(lmin, lmax + 1)), "residual_per_l": per_l, "gain_per_l": gain,
        "alm_in": a_in, "alm_out": a_out,
        "paper": {"relative_residual": PAPER_PIPELINE["residual_fraction"], "map": SOURCE_PLANCK,
                  "source": PAPER_PIPELINE["source"]},
    }


# --- export ---------------------------------------------------------------------------------

def _description(sky: SkyData) -> str:
    name = {SOURCE_WMAP: "WMAP 9-yr ILC", SOURCE_PLANCK: "Planck PR3 Commander"}.get(sky.source_name, sky.source_name)
    return f"{name} at Nside {sky.nside}, monopole and dipole removed"


def sky_map_export(sky: SkyData) -> tuple[dict, dict[str, np.ndarray]]:
    """``sky_map.json`` content and the binary arrays (Float32 LE RING maps, UInt8 mask) by file name.

    When the KQ75 mask was unavailable no mask array is returned and ``mask_file`` is null.
    """
    full = np.asarray(sky.map_uK, dtype="<f4")
    low = np.asarray(sky.lowl_map_uK, dtype="<f4")
    has_mask = sky.mask_fraction is not None

    def _entry(name: str, arr: np.ndarray, desc: str) -> dict:
        r = float(np.percentile(np.abs(arr), 99.5))
        return {"name": name, "file": f"sky_map_{name}.bin", "description": desc,
                "min": float(arr.min()), "max": float(arr.max()), "range": [-r, r]}

    meta = {
        "nside": sky.nside, "ordering": "RING", "frame": "galactic", "unit": "uK",
        "maps": [
            _entry("full", full, _description(sky)),
            _entry("lowl", low, f"{sky.lmin} <= l <= {sky.lmax} reconstruction from the extracted a_lm "
                                f"(beam and pixel window removed)"),
        ],
        "mask_file": "sky_mask.bin" if has_mask else None,
        "mask": {"name": "WMAP KQ75 (9-yr)", "fsky": sky.fsky,
                 "threshold": (sky.meta.get("mask") or {}).get("threshold")} if has_mask else None,
        "source": sky.source_name,
        "computed": True,
    }
    arrays = {"sky_map_full.bin": full, "sky_map_lowl.bin": low}
    if has_mask:
        arrays["sky_mask.bin"] = np.asarray(sky.mask, dtype=np.uint8)
    return meta, arrays


def sky_summary(sky: SkyData) -> dict:
    """JSON-safe metadata of a :class:`SkyData` (no large arrays)."""
    return json.loads(json.dumps({
        "source": sky.source_name, "nside": sky.nside, "lmin": sky.lmin, "lmax": sky.lmax,
        "ell": sky.ell.tolist(), "D_ell": sky.D_ell.tolist(), "fsky": sky.fsky, "meta": sky.meta,
    }, default=float))
