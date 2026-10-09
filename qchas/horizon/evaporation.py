"""Exactly unitary evaporation and Hayden-Preskill recovery, [BH] Sec. 11 and [QG] Theorem R.

The black hole B is n qubits. Before each emission it is scrambled by an exactly
Haar-random unitary; this is applied as a Haar isometry on the Schmidt span of B, which
is equivalent to a Haar unitary on B but costs only O(d_B r^2). One qubit then moves
from B to the radiation R. A diary qubit, maximally entangled with an outside reference,
can be thrown in at any step. Entropies are in bits.
"""

from __future__ import annotations

import numpy as np


def haar_isometry(rng: np.random.Generator, rows: int, cols: int) -> np.ndarray:
    """rows x cols isometry with Haar distribution (QR of a Ginibre matrix, phase-corrected)."""
    z = (rng.standard_normal((rows, cols)) + 1j * rng.standard_normal((rows, cols))) / np.sqrt(2)
    q, r = np.linalg.qr(z)
    d = np.diagonal(r)
    return q * (d / np.abs(d))


def haar_unitary(rng: np.random.Generator, n: int) -> np.ndarray:
    return haar_isometry(rng, n, n)


def _entropy_bits(psi_matrix: np.ndarray) -> float:
    """Entanglement entropy (bits) across the row/column split of a pure-state matrix."""
    m = psi_matrix if psi_matrix.shape[0] <= psi_matrix.shape[1] else psi_matrix.conj().T
    rho = m @ m.conj().T
    w = np.linalg.eigvalsh(rho)
    w = w[w > 1e-14]
    return float(-np.sum(w * np.log2(w)))


def evaporate(n_qubits: int = 20, diary_step: int | None = None, seed: int = 0) -> dict:
    """Emit every qubit of an n-qubit black hole, scrambling exactly before each emission.

    The global state is kept as a tensor psi[b, r, f] over black hole, radiation and
    reference. Returns S(R) after each emission and, if a diary is thrown in before
    emission ``diary_step`` (1-based), the mutual information I(Ref : R) in bits.
    """
    rng = np.random.default_rng(seed)
    psi = np.zeros((2**n_qubits, 1, 1), dtype=complex)
    psi[0, 0, 0] = 1.0
    steps, s_rad, s_bh_cap, i_diary = [], [], [], []
    total = n_qubits + (1 if diary_step is not None else 0)
    for step in range(1, total + 1):
        if diary_step is not None and step == diary_step:
            # diary qubit d enters B, maximally entangled with a fresh reference qubit f
            db, dr, _ = psi.shape
            bell = np.array([[1, 0], [0, 1]], dtype=complex) / np.sqrt(2)  # [d, f]
            psi = np.einsum("brx,df->bdrf", psi[:, :, :1], bell).reshape(db * 2, dr, 2)
        db, dr, dfr = psi.shape
        # scramble B: replace its left singular vectors by a Haar isometry
        mat = psi.reshape(db, dr * dfr)
        u, sv, vh = np.linalg.svd(mat, full_matrices=False)
        keep = sv > 1e-13
        u_new = haar_isometry(rng, db, int(keep.sum()))
        mat = (u_new * sv[keep]) @ vh[keep]
        psi = mat.reshape(db, dr, dfr)
        # emit the last qubit of B into R
        psi = psi.reshape(db // 2, 2, dr, dfr).reshape(db // 2, 2 * dr, dfr)
        db, dr, dfr = psi.shape
        steps.append(step)
        s_rad.append(_entropy_bits(psi.transpose(1, 0, 2).reshape(dr, db * dfr)))
        s_bh_cap.append(float(np.log2(db)))
        if dfr == 2:
            s_ref = _entropy_bits(psi.transpose(2, 0, 1).reshape(dfr, db * dr))
            s_rf = _entropy_bits(psi.transpose(1, 2, 0).reshape(dr * dfr, db))
            i_diary.append(s_ref + s_rad[-1] - s_rf)
        else:
            i_diary.append(0.0)
    s_h = [min(float(s), float(total)) for s in steps]  # Hawking: one bit per emitted qubit
    envelope = [min(a, b) for a, b in zip(s_h, s_bh_cap)]
    return {
        "n_qubits": n_qubits, "diary_step": diary_step, "seed": seed, "steps": steps,
        "s_rad": s_rad, "s_bh": s_bh_cap, "hawking": s_h, "envelope": envelope,
        "page_average": [page_average_entropy(int(round(2 ** sb)), 2 ** (s + 0)) for sb, s in zip(s_bh_cap, steps)],
        "diary_mi": i_diary,
    }


def page_average_entropy(m: int, n: float) -> float:
    """Page's exact average entanglement entropy (bits) of a random pure state on C^m x C^n."""
    m, n = int(min(m, n)), int(max(m, n))
    if m <= 1:
        return 0.0
    k = np.arange(n + 1, m * n + 1)
    return float((np.sum(1.0 / k) - (m - 1) / (2.0 * n)) / np.log(2))


# --- Hayden-Preskill decoding, [QG] Theorem R -------------------------------------------

def _sqrtm_psd(a: np.ndarray) -> np.ndarray:
    w, v = np.linalg.eigh(a)
    return (v * np.sqrt(np.clip(w, 0, None))) @ v.conj().T


def _fidelity(rho: np.ndarray, sigma: np.ndarray) -> float:
    """Root fidelity ||sqrt(rho) sqrt(sigma)||_1."""
    sv = np.linalg.svd(_sqrtm_psd(rho) @ _sqrtm_psd(sigma), compute_uv=False)
    return float(np.sum(sv))


def decoding_infidelity(n_bh: int = 7, k: int = 2, s_values=range(0, 10), samples: int = 6,
                        old: bool = True, dynamics: str = "haar", t_eth: float = 10.0,
                        seed: int = 0) -> dict:
    """Upper bound 1 - F(rho_MB', tau_M x rho_B')^2 on the decoding infidelity.

    An old black hole (``old=True``) is n_bh qubits maximally entangled with n_bh qubits of
    early radiation E; a young one starts pure. A k-qubit diary D, maximally entangled
    with reference M, is thrown in; B+D is scrambled (Haar unitary, or e^{-iHt} with H from
    the GUE rescaled to [-1, 1]); then s qubits are emitted. B' is what remains.
    """
    rng = np.random.default_rng(seed)
    nb = n_bh + k
    out = {s: [] for s in s_values}
    for _ in range(samples):
        db, dd, dm = 2**n_bh, 2**k, 2**k
        de = 2**n_bh if old else 1
        # psi[b, e] for the black hole and early radiation
        bh = np.eye(db, dtype=complex)[:, :de] / np.sqrt(de) if old else np.eye(db, 1, dtype=complex)
        diary = np.eye(dd, dtype=complex) / np.sqrt(dd)  # [d, m]
        psi = np.einsum("be,dm->bdem", bh, diary).reshape(db * dd, de, dm)
        if dynamics == "haar":
            u = haar_unitary(rng, 2**nb)
        else:
            g = (rng.standard_normal((2**nb,) * 2) + 1j * rng.standard_normal((2**nb,) * 2)) / 2
            h = g + g.conj().T
            w, v = np.linalg.eigh(h)
            w = 2 * (w - w.min()) / (w.max() - w.min()) - 1
            u = (v * np.exp(-1j * w * t_eth)) @ v.conj().T
        psi = np.einsum("ij,jem->iem", u, psi)
        for s in s_values:
            d_rest = 2 ** (nb - s)
            # B = B' (first nb - s qubits) x emitted (last s qubits)
            t = psi.reshape(d_rest, 2**s, de, dm)
            rho_mb = np.einsum("aiem,bien->manb", t, t.conj()).reshape(dm * d_rest, dm * d_rest)
            rho_b = np.einsum("aiem,biem->ab", t, t.conj())
            target = np.kron(np.eye(dm) / dm, rho_b)
            f = _fidelity(rho_mb, target)
            out[s].append(1.0 - f**2)
    return {"n_bh": n_bh, "k": k, "old": old, "dynamics": dynamics, "t": t_eth,
            "s": list(s_values), "infidelity": [float(np.mean(out[s])) for s in s_values],
            "bound": [float(min(1.0, 2.0 ** (1 + k - s))) for s in s_values]}
