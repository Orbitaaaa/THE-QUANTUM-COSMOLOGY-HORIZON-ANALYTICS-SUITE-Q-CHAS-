"""Round 11 (G11-12): greybody tables for model U (Hayward, D = 10) and Tangherlini (n = 6) -> DATA/round11_greybody.json."""
import json
import sys
import time

from _common import DATA
from babyuniverse import greybody11 as G

if __name__ == "__main__":
    nw = int(sys.argv[sys.argv.index("--workers") + 1]) if "--workers" in sys.argv else 4
    t = time.time()
    out = dict(hayward=G.emission_tables(10, (1.05, 1.5, 3.0, 10.0), "hayward", workers=nw))
    print("hayward", time.time() - t, flush=True)
    out["gr"] = G.emission_tables(10, (1.0,), "gr", workers=nw)
    print("gr", time.time() - t, flush=True)
    json.dump(out, open(f"{DATA}/round11_greybody.json", "w"))
    print("wrote DATA/round11_greybody.json", flush=True)
