"""Collect every paper number Q-CHAS reproduces into one table (``validation.json``).

Reads the bundles already written by :mod:`qchas.export` and lists, per quantity, our value,
the printed paper value, the relative difference and a verdict. Verdicts are honest:
``match`` within the stated tolerance, ``differs`` otherwise, with the diagnosis that the
analysis module recorded (for example a paper typo or a grid effect).
"""

from __future__ import annotations

import json
from pathlib import Path


def _load(out_dir: Path, name: str):
    p = out_dir / name
    return json.loads(p.read_text()) if p.exists() else None


class _Rows:
    def __init__(self):
        self.rows: list[dict] = []

    def add(self, group: str, quantity: str, ours, paper, source: str, tol: float = 0.02, note: str = ""):
        if ours is None or paper is None:
            return
        try:
            o, p = float(ours), float(paper)
        except (TypeError, ValueError):
            return
        rel = abs(o - p) / abs(p) if p != 0 else abs(o - p)
        # agreement to the printed precision also counts as a match
        printed = _printed_match(o, paper)
        verdict = "match" if (rel <= tol or printed) else "differs"
        self.rows.append({"group": group, "quantity": quantity, "ours": o, "paper": p,
                          "rel_diff": rel, "verdict": verdict, "source": source, "note": note})


def _printed_match(ours: float, paper) -> bool:
    """True when ``ours`` rounds to the paper value at the precision the paper printed."""
    s = repr(paper) if not isinstance(paper, str) else paper
    try:
        p = float(s)
    except ValueError:
        return False
    if "e" in s.lower():
        mant = s.lower().split("e")[0]
        digits = len(mant.split(".")[1]) if "." in mant else 0
        exp = int(s.lower().split("e")[1])
        step = 10.0 ** (exp - digits)
    else:
        digits = len(s.split(".")[1]) if "." in s else 0
        step = 10.0 ** (-digits)
    return abs(ours - p) <= 0.5 * step * (1 + 1e-9)


def _group_qes(r: _Rows, horizon: dict):
    for row in horizon["qes"]["rows"]:
        tag = f"M={row['mass']:g}, c={row['c']:g}, r_b={row['rb_over_rs']:g} r_s"
        p = row["paper"]
        r.add("QES ([BH] Sec. 9)", f"depth (r_s - r_QES)/r_s, {tag}", row["depth"], p["depth"], "[BH] Sec. 9", tol=0.005)
        r.add("QES ([BH] Sec. 9)", f"lag u_b - v_QES, {tag}", row["lag"], p["lag"], "[BH] Sec. 9", tol=0.001)
        r.add("QES ([BH] Sec. 9)", f"S_gen - S_BH, {tag}", row["s_gen_minus_s_bh"], p["s_gen_minus_s_bh"], "[BH] Sec. 9", tol=0.01)


def _group_atmosphere(r: _Rows, horizon: dict):
    a = horizon["atmosphere"]
    r.add("Quantum atmosphere ([BH] Sec. 10)", "area-law coefficient S/(n+1/2)^2", a["area_coeff"], 0.295, "[BH] Sec. 10", tol=0.005)
    for g, ours, paper in zip(a["gaps"], a["share"], a["paper_share"]):
        if g and paper is not None:
            r.add("Quantum atmosphere ([BH] Sec. 10)", f"share beyond gap g={g}", ours, paper, "[BH] Sec. 10", tol=0.01)


def _group_page(r: _Rows, horizon: dict):
    pg = horizon["page"]
    r.add("Page transfer law ([BH] Fig. 1)", "Page time t_P / t_evap", pg["page_time"], 0.476, "[BH] Fig. 1", tol=0.002)
    r.add("Page transfer law ([BH] Fig. 1)", "mass fraction at the Page time", pg["page_mass_fraction"], 0.806, "[BH] Fig. 1", tol=0.002)
    r.add("Crossed product ([QG] Theorem U)", "area gap a0 / l_P^2", horizon["capacity"]["min_area_lp2"], 5.97, "[QG] Theorem U", tol=0.002)


def _group_evaporation(r: _Rows, horizon: dict):
    ev = horizon["evaporation"]
    paper = ev.get("paper", {})
    early, late = ev.get("early_diary", {}), ev.get("late_diary", {})
    if early.get("n_qubits") != paper.get("n_qubits"):
        return
    g = "Unitary evaporation ([BH] Sec. 11)"
    r.add(g, "maximum radiation entropy (bits)", early.get("max_bits"), paper.get("max_bits"), "[BH] Sec. 11", tol=0.002)
    r.add(g, "step of the maximum", early.get("max_step"), paper.get("max_step"), "[BH] Sec. 11", tol=0)
    r.add(g, "largest deviation from min{S_H, S_BH} (bits)", early.get("max_dev"), paper.get("max_dev"), "[BH] Sec. 11", tol=0.01)
    r.add(g, "early diary (in at 5) recovered at step", early.get("recovered_step"), paper.get("early_recovered_step"), "[BH] Sec. 11", tol=0)
    r.add(g, "late diary (in at 14) recovered at step", late.get("recovered_step"), paper.get("late_recovered_step"), "[BH] Sec. 11", tol=0)


def _group_decoding(r: _Rows, horizon: dict):
    dec = horizon["decoding"]
    paper = dec.get("paper", {})
    haar = dec.get("haar_old", {})
    p_haar = paper.get("haar") or paper.get("haar_old")
    if not p_haar:
        return
    for s, ours, p in zip(haar.get("s", []), haar.get("infidelity", []), p_haar):
        if p is not None and p > 1e-10:
            r.add("Hayden-Preskill decoding ([QG] Theorem R)", f"Haar infidelity bound after s={s} emissions", ours, p,
                  "[QG] Sec. 15 table", tol=0.25, note="mean of 6 Haar samples; sampling noise ~10-20 %")


def _group_lattice(r: _Rows, lattice: dict):
    g = "Casimir lattice ([TOP] Prop. E.2)"
    for s, c in lattice["critical"].items():
        r.add(g, f"Z(2) at {s}", c["Z"], c["Z_paper"], "[TOP] Prop. E.2", tol=1e-12)
        seen = set()
        for ours, paper in zip(sorted(c["hessian"]), sorted(c["hessian_paper"])):
            if round(paper, 4) in seen:
                continue
            seen.add(round(paper, 4))
            note = ("printed value violates the exact identity tr H = 2Z/3; analytic and finite-difference "
                    "Hessians agree on 3.39504 (likely typo)") if (s == "BCC" and abs(paper - 3.3953) < 1e-6) else ""
            r.add(g, f"Hessian eigenvalue at {s}", ours, paper, "[TOP] Prop. E.2", tol=5e-4, note=note)
    fl = lattice.get("flows", {})
    if fl.get("fcc_fraction") is not None:
        r.add("Shape flow ([TOP] Sec. 6.3)", "fraction of flows ending at FCC", fl["fcc_fraction"], 0.79,
              "[TOP] Prop. H.3 / Sec. 14", tol=0.1, note=f"N = {fl.get('n')}, kappa = {fl.get('kappa')}; paper N = 240, kappa unstated")


def _group_circles(r: _Rows, topo: dict):
    want = {("SC", 0.985): 3, ("BCC", 0.985): 4, ("FCC", 0.985): 6, ("SC", 0.8): 3, ("BCC", 0.8): 7, ("FCC", 0.8): 6}
    for t in topo.get("circles", {}).get("tables", []):
        key = (t["shape"], round(t["d_over_D"], 3))
        if key in want:
            r.add("Circles in the sky ([TOP] Sec. 7)", f"circle pairs, {key[0]} at d/D = {key[1]}", t["n_pairs"], want[key],
                  "[TOP] Sec. 7 table, Fig. 4", tol=0)
            if key[1] == 0.985:
                r.add("Circles in the sky ([TOP] Sec. 7)", f"circle radius (deg), {key[0]} at d/D = 0.985", t["radius_deg"], 9.9,
                      "[TOP] Sec. 7 table", tol=0.01)


def _group_kl(r: _Rows, cmb: dict):
    kl = cmb.get("kl")
    if not kl:
        return
    for shape in ("cubic", "BCC", "FCC"):
        for d, ours in zip(kl["d_over_D"], kl[shape]):
            p = kl["paper"][shape].get(f"{d:g}") or kl["paper"][shape].get(str(d))
            r.add("Detectability KL ([TOP] Theorem P)", f"KL(l<=30), {shape} at d/D = {d:g}", ours, p, "[TOP] Theorem P table", tol=0.02)
    sec8 = cmb.get("sec8")
    if sec8:
        r.add("Large-angle CMB ([TOP] Sec. 8)", "LCDM D_2 (uK^2)", sec8["lcdm_D"][0], sec8.get("lcdm_D2_paper", 1016), "[TOP] Sec. 8", tol=0.005)
        for row in sec8.get("rows", []):
            tag = f"{row['shape']} at d/D = {row['d_over_D']:g}"
            p = row.get("paper", {})
            r.add("Large-angle CMB ([TOP] Sec. 8)", f"D_2 (uK^2), {tag}", row.get("D2"), p.get("D2"), "[TOP] Sec. 8 table", tol=0.005)
            r.add("Large-angle CMB ([TOP] Sec. 8)", f"KL(l<=30), {tag}", row.get("KL") or row.get("kl"), p.get("KL") or p.get("kl"),
                  "[TOP] Sec. 8 table", tol=0.02)


def _group_resonance(r: _Rows, res: dict):
    for row in res.get("kl_table", []):
        tag = f"cubic d/D = {row['d_over_D']:g}, shell |m|^2 = {row['shell']} ({row['n_vectors']} vectors), f = {row['f']:g}"
        r.add("Topological resonance ([TOP] Sec. 11.6)", f"KL with line, {tag}", row["kl"], row["paper"]["kl"], "[TOP] Sec. 11.6 table", tol=0.02)
    seen = set()
    for row in res.get("kl_table", []):
        if row["d_over_D"] in seen:
            continue
        seen.add(row["d_over_D"])
        r.add("Topological resonance ([TOP] Sec. 11.6)", f"information gain at d/D = {row['d_over_D']:g}", row["gain"],
              row["paper"]["gain"], "[TOP] Sec. 11.6 table", tol=0.08, note="paper prints the gain rounded to an integer")


def _group_cmb_stats(r: _Rows, cmb: dict):
    st = cmb.get("stats", {})
    pc = st.get("paper_check") or {}
    g = "Low-l anomalies ([TOP] Sec. 8)"
    if pc:
        r.add(g, "LCDM median S_1/2 (uK^4)", pc.get("S12_lcdm_median"), 34576, "[TOP] Sec. 8 table", tol=0.03,
              note="Monte Carlo median; seed scatter ~1 %")
        r.add(g, "P(S_1/2 <= 6694 uK^4) under LCDM", pc.get("p_S12_planck_value"), 0.076, "[TOP] Sec. 8 table", tol=0.06)
        r.add(g, "P(C_2 <= Planck's 226 uK^2) under LCDM", pc.get("p_C2_planck_value"), 0.047, "[TOP] Sec. 8 table", tol=0.06)


def _group_bounce(r: _Rows, bounce: dict):
    g = "Quantum Big Bang ([QG] Theorem Lambda)"
    for run in bounce.get("runs", []):
        p, fd = run.get("paper"), run.get("paper_grid")
        if not p or "peak_a_inv2" not in p or not fd:
            continue
        a = run["a_max"]
        r.add(g, f"peak <a^-2>, a_max = {a:g} (paper's finite-difference grid)", fd.get("peak_a_inv2"), p["peak_a_inv2"],
              "[QG] Theorem Lambda table", tol=0.005)
        r.add(g, f"<a>_min, a_max = {a:g} (paper's finite-difference grid)", fd.get("a_min"), p["a_min"],
              "[QG] Theorem Lambda table", tol=0.002)
        r.add(g, f"FWHM DeltaT, a_max = {a:g} (paper's finite-difference grid)", fd.get("fwhm_T"), p["fwhm_T"],
              "[QG] Theorem Lambda table", tol=0.01)
        r.add(g, f"peak <a^-2>, a_max = {a:g} (exact solution)", run.get("peak_a_inv2"), p["peak_a_inv2"],
              "[QG] Theorem Lambda table", tol=0.005,
              note="exact closed-form solution; the printed table carries an O(h) grid bias at a = 0")


def _group_gw(r: _Rows, gw: dict):
    g = "Kerr ringdown and echoes ([QG] Sec. 19)"
    for row in gw.get("paper_events", []):
        r.add(g, f"f_220 (Hz), {row['event']}", row.get("f220_ours"), row.get("f220"), "[QG] Sec. 19 table", tol=0.01)
        r.add(g, f"echo delay (ms), {row['event']}", row.get("dt_ms_ours"), row.get("dt_ms"), "[QG] Sec. 19 table", tol=0.012)
    chk = gw.get("echo_bound_check", {})
    if chk.get("log10_R2_ours") is not None:
        paper = chk.get("log10_R2") or chk.get("log10_R2_paper")
        r.add(g, "echo power bound log10|R|^2 for a 62 Msun remnant", chk["log10_R2_ours"], paper, "[QG] Sec. 18", tol=0.005)


def collect(out_dir: Path) -> dict:
    out_dir = Path(out_dir)
    r = _Rows()
    bundles = {n: _load(out_dir, f"{n}.json") for n in ("horizon", "lattice", "topology", "cmb", "resonance", "bounce", "gw")}
    groups = [
        (_group_qes, "horizon"), (_group_atmosphere, "horizon"), (_group_page, "horizon"),
        (_group_evaporation, "horizon"), (_group_decoding, "horizon"), (_group_lattice, "lattice"),
        (_group_circles, "topology"), (_group_kl, "cmb"), (_group_cmb_stats, "cmb"),
        (_group_resonance, "resonance"), (_group_bounce, "bounce"), (_group_gw, "gw"),
    ]
    skipped = []
    for fn, name in groups:
        if bundles.get(name) is None:
            skipped.append(f"{fn.__name__[7:]}: {name}.json missing")
            continue
        try:
            fn(r, bundles[name])
        except (KeyError, TypeError, IndexError) as exc:
            skipped.append(f"{fn.__name__[7:]}: {type(exc).__name__}: {exc}")
    n_match = sum(1 for x in r.rows if x["verdict"] == "match")
    return {
        "rows": r.rows,
        "summary": {"total": len(r.rows), "match": n_match, "differs": len(r.rows) - n_match},
        "skipped": skipped,
        "note": ("Every row compares a number Q-CHAS computes from scratch with the value printed in the "
                 "source paper. 'match' means agreement within the stated tolerance or to the printed precision."),
        "computed": True,
    }


__all__ = ["collect"]
