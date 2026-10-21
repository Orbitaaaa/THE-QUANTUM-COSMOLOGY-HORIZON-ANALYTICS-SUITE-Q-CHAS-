"""Run every symbolic derivation and write DATA/derivation_log.txt."""
import os
import time

from _common import DATA, ROOT  # noqa: F401
from babyuniverse import symbolic

lines = ["Symbolic derivation log (sympy %s)" % symbolic.sp.__version__, ""]
t0 = time.time()
res = symbolic.run_all()
for name, d in res.items():
    lines.append(f"[verified] {name}")
    for k, v in d.items():
        lines.append(f"    {k} = {v}")
    lines.append("")
lines.append(f"total time {time.time() - t0:.1f} s")
path = os.path.join(DATA, "derivation_log.txt")
open(path, "w").write("\n".join(lines))
print("wrote", os.path.relpath(path, ROOT))
