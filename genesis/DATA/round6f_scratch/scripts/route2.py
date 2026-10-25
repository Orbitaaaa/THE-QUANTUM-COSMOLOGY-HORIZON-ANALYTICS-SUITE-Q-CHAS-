import numpy as np, warnings
warnings.simplefilter("ignore")
from babyuniverse import condensate as C
for lf in (1e-6, 3e-6):
    o = C.route_condensate(lam_f=lf)
    res = o["res"]; R, P, Cc, end = res["r"], res["phi"], o["cls"], res["end"]
    h = 0.04
    for i in (10, 25, 50, 80, 150, 300, 600, 749):
        js = np.where(np.isfinite(R[i]))[0]; j = js.max()
        tv_norm_far = np.sum((P[i] > 0) & (Cc[i] >= 0) & (R[i] > 5 * o["r_N"]))
        print(f"lam={lf} u={i*h:5.2f} end_code={end[i]} v_end={j*h:7.2f} r_end={R[i,j]:9.3f} phi_end={P[i,j]:+.4f} cls_end={Cc[i,j-1]} "
              f"TV_normal_far_pts={tv_norm_far} anti_FV_big={np.sum((Cc[i]==1)&(P[i]<0)&(R[i]>5*o['r_N']))}")
