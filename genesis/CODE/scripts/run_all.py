"""Reproduce every data file and figure (Deliverable 1, Phases II-VI, roadmap, round 3)."""
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = ["derivation_log.py", "branch_atlas.py", "tunnelling_analysis.py", "singularity_audit.py",
           "thickwall_study.py", "quantum_checks.py", "physical_scales.py",
           "phase2_verification.py", "phase2_regimeB_scan.py", "phase3_analysis.py", "phase4_dynamics.py", "phase5_6_analysis.py", "roadmap_analysis.py", "sds_atlas_figure.py",
           "round3_analysis.py", "round3_stage1.py", "round6_atlas.py", "round5_7.py", "round6b_gsl.py", "round6c_throat.py", "round6c_nariai.py",
           "round6c_answers.py", "round6d.py", "round6e.py", "round6f.py", "round6g.py", "round6i.py", "round6j.py", "round6k.py", "round6l.py", "round6m.py", "round6n.py", "round6o.py", "round6p.py", "round7.py", "round8.py", "round9.py", "round10_open.py", "round10a.py", "round10.py", "machine10.py"]
for s in SCRIPTS:
    t0 = time.time()
    print(f"=== {s}", flush=True)
    r = subprocess.run([sys.executable, os.path.join(HERE, s)], cwd=HERE)
    if r.returncode != 0:
        sys.exit(f"{s} failed")
    print(f"=== {s} done in {time.time() - t0:.0f} s", flush=True)
