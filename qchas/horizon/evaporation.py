"""Exactly unitary evaporation and Hayden-Preskill recovery, [BH] Sec. 11 and [QG] Theorem R.

The black hole B is n qubits in a product state. Before each emission it is scrambled by
an exactly Haar-random unitary, applied as a Haar isometry on the Schmidt span of B (the
global state is a d_B x d_R d_Ref matrix; replacing its left singular vectors by a Haar
isometry is exactly a Haar unitary on B, [BH] Sec. 17). One qubit then moves from B to
the radiation R. A diary qubit, maximally entangled with an outside reference, can be
thrown into B before any emission. Entropies are in bits.

Theorem R is tested as in [QG] Sec. 15: an old black hole of n_B qubits maximally
entangled with early radiation E swallows a k-qubit diary D (entangled with M), B+D is
scrambled by a Haar unitary or by e^{-iHt} with H from the GUE rescaled to [-1, 1], and
s qubits of the original black hole are emitted first (the diary's own qubits last).
"""

from __future__ import annotations

import math

import numpy as np
from scipy.special import digamma

#: [BH] Sec. 11 (20 qubits, one diary qubit)
PAPER_EVAPORATION = {
    "n_qubits": 20, "max_bits": 10.28, "max_step": 11, "max_dev": 0.72, "final_bits_with_diary": 1.0,
    "early_diary_step": 5, "early_recovered_step": 13, "late_diary_step": 14, "late_recovered_step": 17,
    "recovery_threshold_bits": 1.9, "source": "[BH] Sec. 11",
}

#: [QG] Theorem R table (p. 28): mean of 6 samples, n_B = 7, k = 2; None where printed as ~1e-16
PAPER_DECODING = {
    "n_bh": 7, "k": 2, "samples": 6, "s": list(range(10)),
    "bound": [1, 1, 1, 1, 0.5, 0.25, 0.125, 0.0625, 0.0312, 0.0156],
    "haar_old": [0.937, 0.762, 0.271, 0.0601, 0.0148, 0.00366, 0.000903, 0.000209, 4.13e-5, 5.55e-16],
    "eth": {
        "1": [0.937, 0.841, 0.648, 0.571, 0.554, 0.550, 0.549, 0.549, 0.336, 5.55e-16],
        "3": [0.937, 0.763, 0.274, 0.0656, 0.0208, 0.00976, 0.00697, 0.00629, 0.00153, 3.33e-16],
        "10": [0.937, 0.762, 0.270, 0.0601, 0.0146, 0.00364, 0.000913, 0.000208, 4.19e-5, 7.77e-16],
        "100": [0.937, 0.762, 0.270, 0.0599, 0.0148, 0.00370, 0.000907, 0.000220, 4.42e-5, 3.33e-16],
    },
    "haar_young": [0.937, 0.936, 0.931, 0.908, 0.825, 0.530, 0.120, 0.0288, 0.00580, 1.11e-16],
    "source": "[QG] Theorem R table (Sec. 15)",
}


# --- random unitaries ------------------------------------------------------------------

def haar_isometry(rng: np.random.Generator, rows: int, cols: int) -> np.ndarray:
    """rows x cols isometry with Haar distribution (QR of a Ginibre matrix, phase-corrected)."""
    z = (rng.standard_normal((rows, cols)) + 1j * rng.standard_normal((rows, cols))) / np.sqrt(2)
    q, r = np.linalg.qr(z)
    d = np.diagonal(r)
    return q * (d / np.abs(d))


def haar_unitary(rng: np.random.Generator, n: int) -> np.ndarray:
    return haar_isometry(rng, n, n)


def gue_unitary(rng: np.random.Generator, n: int, t: float) -> np.ndarray:
    """U = exp(-iHt), H from the GUE with its spectrum mapped affinely onto [-1, 1]."""
    g = (rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))) / 2
    w, v = np.linalg.eigh(g + g.conj().T)
    w = 2.0 * (w - w[0]) / (w[-1] - w[0]) - 1.0
    return (v * np.exp(-1j * w * t)) @ v.conj().T


# --- entropies -------------------------------------------------------------------------

def _entropy_bits(psi_matrix: np.ndarray) -> float:
    """Entanglement entropy (bits) across the row/column split, from the smaller Gram matrix."""
    m = psi_matrix if psi_matrix.shape[0] <= psi_matrix.shape[1] else psi_matrix.T
    w = np.linalg.eigvalsh(m @ m.conj().T)
    w = w[w > 1e-14]
    return float(-np.sum(w * np.log2(w)))


def page_average_entropy(m: float, n: float) -> float:
    """Page's exact average entropy (bits) of the smaller factor of a random pure state on C^m x C^n.

    S_{m,n} = sum_{k=n+1}^{mn} 1/k - (m - 1)/(2n) nats for m <= n (Page 1993), with the
    harmonic sum written as digamma(mn + 1) - digamma(n + 1).
    """
    m, n = float(min(m, n)), float(max(m, n))
    if m <= 1.0:
        return 0.0
    return float((digamma(m * n + 1.0) - digamma(n + 1.0) - (m - 1.0) / (2.0 * n)) / math.log(2.0))


# --- exactly unitary evaporation, [BH] Sec. 11 ------------------------------------------

def _scramble(rng: np.random.Generator, mat: np.ndarray) -> np.ndarray:
    """Apply a Haar unitary U_B to the rows of the d_B x c state matrix without forming U_B.

    For d_B <= c this is U_B mat. Otherwise mat = Q R (thin QR) and U_B Q is a Haar
    isometry, so U_B mat has the law of W R with W a fresh d_B x c Haar isometry.
    """
    db, c = mat.shape
    if db <= c:
        return haar_unitary(rng, db) @ mat
    r = np.linalg.qr(mat, mode="r")
    return haar_isometry(rng, db, c) @ r


def evaporate(n_qubits: int = 20, diary_step: int | None = None, seed: int = 0,
              recovery_bits: float = PAPER_EVAPORATION["recovery_threshold_bits"]) -> dict:
    """Emit every qubit of an n-qubit black hole, scrambling exactly before each emission.

    The global state is a tensor psi[b, r, f] over black hole, radiation and reference. If
    ``diary_step`` is given, the diary is thrown into B after that many emissions (the
    convention of [BH] Figure 2: "thrown in at step 14" leaves I(Ref : R) = 0 at 14), so the
    run has n + 1 emissions. Per step it records S(R), the Hawking entropy S_H (one bit per
    emitted qubit), the capacity log2 d_B, the reference entropy, the envelope
    min{S_H, log2 d_B + S(Ref)} of Theorem E(i) (S(R) = S(B Ref)), Page's exact average for
    the same dimensions, and the diary's mutual information I(Ref : R).
    """
    rng = np.random.default_rng(seed)
    psi = np.zeros((2**n_qubits, 1, 1), dtype=complex)
    psi[0, 0, 0] = 1.0
    total = n_qubits + (1 if diary_step is not None else 0)
    rec: dict[str, list] = {k: [] for k in ("steps", "s_rad", "hawking", "s_bh", "s_ref", "envelope",
                                            "page_average", "diary_mi")}
    for step in range(1, total + 1):
        if diary_step is not None and step == diary_step + 1:
            db, dr, _ = psi.shape
            bell = np.eye(2, dtype=complex) / np.sqrt(2)  # [d, f]
            psi = np.einsum("brx,df->bdrf", psi[:, :, :1], bell).reshape(db * 2, dr, 2)
        db, dr, dfr = psi.shape
        psi = _scramble(rng, psi.reshape(db, dr * dfr)).reshape(db // 2, 2 * dr, dfr)  # last qubit of B -> R
        db, dr, dfr = psi.shape
        s_r = _entropy_bits(psi.transpose(1, 0, 2).reshape(dr, db * dfr))
        s_ref = _entropy_bits(psi.transpose(2, 0, 1).reshape(dfr, db * dr)) if dfr > 1 else 0.0
        if dfr > 1:
            s_rf = _entropy_bits(psi.transpose(1, 2, 0).reshape(dr * dfr, db))
            mi = s_ref + s_r - s_rf
        else:
            mi = 0.0
        s_bh = math.log2(db)
        n_ref = math.log2(dfr)
        rec["steps"].append(step)
        rec["s_rad"].append(s_r)
        rec["hawking"].append(float(step))
        rec["s_bh"].append(s_bh)
        rec["s_ref"].append(s_ref)
        rec["envelope"].append(min(float(step), s_bh + n_ref))
        rec["page_average"].append(page_average_entropy(2.0**step, 2.0 ** (s_bh + n_ref)))
        rec["diary_mi"].append(mi)
    dev = np.array(rec["envelope"]) - np.array(rec["s_rad"])
    k_max = int(np.argmax(rec["s_rad"]))
    recovered = None
    if diary_step is not None:
        recovered = next((s for s, i in zip(rec["steps"], rec["diary_mi"]) if i >= recovery_bits), None)
    return {
        "n_qubits": n_qubits, "diary_step": diary_step, "seed": seed, **rec,
        "max_bits": rec["s_rad"][k_max], "max_step": rec["steps"][k_max],
        "max_dev": float(dev.max()), "max_dev_step": rec["steps"][int(np.argmax(dev))],
        "final_bits": rec["s_rad"][-1], "recovered_step": recovered,
        "recovery_threshold_bits": recovery_bits, "computed": True,
    }


# --- Hayden-Preskill decoding, [QG] Theorem R -------------------------------------------

def _infidelity_after(psi: np.ndarray, n_total: int, s: int) -> float:
    """1 - F(rho_MB', tau_M x rho_B')^2 for psi[i, e, m] with i = (b', emitted s qubits).

    With sigma = tau_M x rho_B', F = ||(1 x sqrt(rho_B')) X||_1 / sqrt(d_M), where X is the
    (M B') x (R' E) matrix of the pure state; the trace norm comes from the smaller Gram matrix.
    """
    _, d_e, d_m = psi.shape
    d_bp = 2 ** (n_total - s)
    x = psi.reshape(d_bp, 2**s, d_e, d_m).transpose(3, 0, 1, 2).reshape(d_m, d_bp, 2**s * d_e)
    xb = x.transpose(1, 0, 2).reshape(d_bp, -1)
    rho_b = xb @ xb.conj().T
    w, v = np.linalg.eigh(rho_b)
    sqrt_rho = (v * np.sqrt(np.clip(w, 0.0, None))) @ v.conj().T
    y = np.matmul(sqrt_rho[None], x).reshape(d_m * d_bp, -1)
    gram = y.conj().T @ y if y.shape[1] <= y.shape[0] else y @ y.conj().T
    f = float(np.sum(np.sqrt(np.clip(np.linalg.eigvalsh(gram), 0.0, None)))) / math.sqrt(d_m)
    return 1.0 - f * f


def decoding_infidelity(n_bh: int = 7, k: int = 2, s_values=None, samples: int = 6,
                        old: bool = True, dynamics: str = "haar", t_eth: float = 10.0,
                        seed: int = 0) -> dict:
    """Mean and spread of 1 - F(rho_MB', tau_M x rho_B')^2, the upper bound on 1 - F_dec.

    B+D is indexed (d, b) with the diary most significant, and emission takes the least
    significant qubits, so the original black-hole qubits leave first. An old black hole
    is maximally entangled with n_bh qubits of early radiation E; a young one starts in
    |0...0>. ``dynamics`` is "haar" or "eth" (U = e^{-iHt}, GUE spectrum on [-1, 1]).
    """
    rng = np.random.default_rng(seed)
    n_total = n_bh + k
    d_b, d_d, dim = 2**n_bh, 2**k, 2**(n_bh + k)
    s_values = list(range(n_total + 1) if s_values is None else s_values)
    if max(s_values) > n_total:
        raise ValueError(f"cannot emit more than the {n_total} qubits of B+D")
    vals = np.empty((samples, len(s_values)))
    for i in range(samples):
        u = haar_unitary(rng, dim) if dynamics == "haar" else gue_unitary(rng, dim, t_eth)
        if old:
            # |Phi_BE>|Phi_DM>: column (m, e) of U, normalised
            psi = u.reshape(dim, d_d, d_b).transpose(0, 2, 1) / math.sqrt(d_b * d_d)
        else:
            psi = u[:, ::d_b].reshape(dim, 1, d_d) / math.sqrt(d_d)
        for j, s in enumerate(s_values):
            vals[i, j] = _infidelity_after(psi, n_total, s)
    return {"n_bh": n_bh, "k": k, "old": old, "dynamics": dynamics, "t": t_eth if dynamics == "eth" else None,
            "samples": samples, "seed": seed, "s": s_values,
            "infidelity": vals.mean(axis=0).tolist(), "std": vals.std(axis=0, ddof=1).tolist() if samples > 1 else None,
            "bound": [float(min(1.0, 2.0 ** (1 + k - s))) for s in s_values], "computed": True}


def decoding_table(samples: int = 6, seed: int = 0, n_bh: int = 7, k: int = 2) -> dict:
    """Every column of the [QG] Theorem R table, with independent seeds per column."""
    out = {
        "haar_old": decoding_infidelity(n_bh, k, samples=samples, seed=seed),
        "eth": {str(t): decoding_infidelity(n_bh, k, samples=samples, dynamics="eth", t_eth=float(t),
                                            seed=seed + 10 + i)
                for i, t in enumerate((1, 3, 10, 100))},
        "haar_young": decoding_infidelity(n_bh, k, samples=samples, old=False, seed=seed + 1),
        "paper": PAPER_DECODING,
    }
    return out


def to_json(n_qubits: int = 20, samples: int = 6, seed: int = 0) -> dict:
    """The "evaporation" and "decoding" parts of horizon.json."""
    early, late = PAPER_EVAPORATION["early_diary_step"], PAPER_EVAPORATION["late_diary_step"]
    if n_qubits != PAPER_EVAPORATION["n_qubits"]:
        early, late = max(1, round(n_qubits / 4)), round(0.7 * n_qubits)
    return {
        "evaporation": {
            "no_diary": evaporate(n_qubits, None, seed),
            "early_diary": evaporate(n_qubits, early, seed + 1),
            "late_diary": evaporate(n_qubits, late, seed + 2),
            "paper": PAPER_EVAPORATION,
        },
        "decoding": decoding_table(samples, seed),
        "claims": CLAIMS,
    }


CLAIMS = [
    {"statement": "S(R) <= min{S_H, A/4G}; equality up to e^{-|S_H - A/4G|} for 2-design scrambling",
     "status": "Theorem, given capacity and H*", "source": "[BH] Theorem E"},
    {"statement": "Exactly unitary 20-qubit evaporation peaks at 10.28 bits at the Page time and returns 1 bit (the diary's entanglement) at the end",
     "status": "Computation", "source": "[BH] Sec. 11"},
    {"statement": "A diary thrown in after the Page time is decoded 3 emissions later (Hayden-Preskill)",
     "status": "Computation", "source": "[BH] Sec. 11, Theorem E(iv)"},
    {"statement": "E[1 - F_dec] <= 2^{1+k-s} for 2-design dynamics; GUE e^{-iHt} matches Haar from t = 10",
     "status": "Theorem for 2-design scrambling; tested with the ETH dynamics", "source": "[QG] Theorem R"},
]
