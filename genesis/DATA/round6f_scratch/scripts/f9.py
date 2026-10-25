import numpy as np, warnings
warnings.simplefilter("ignore")
from scipy.integrate import solve_ivp
from babyuniverse import throat_evolution as TE, condensate as C
rho_f = 3/900/(8*np.pi); v = 0.0015; lam = 4*rho_f/(4*v**4)
pot = TE.Potential(lam, v, rho_f)
dU = lambda p: lam*p*(p*p - v*v) - rho_f*(0.75 - 0.75*(p/v)**2)/v
sol = solve_ivp(lambda t, y: [y[1], -dU(y[0])], (0, 4), [0.5*v, 0.0], max_step=0.0005, rtol=1e-11, atol=1e-16, dense_output=True)
ph = lambda t: sol.sol(np.asarray(t, float))[0]
h = 0.01; N = 400
a, b, _ = TE.characteristic_data(0.0, 0.0, 60.0, -0.5, 0.5, h, N, N, phi_v=lambda x: ph(np.asarray(x)/2), phi_u=lambda x: ph(np.asarray(x)/2), G=1e-10)
r1 = TE.evolve(a, b, 0.0, 0.0, pot, h, G=1e-10)
z = lambda x: 0*np.asarray(x, float)
A, B, _ = C.condensate_data(0.0, 0.0, 60.0, -0.5, 0.5, h, N, N, z, z, phi_v=lambda x: ph(np.asarray(x)/2), phi_u=lambda x: ph(np.asarray(x)/2), G=1e-10, frame="corner")
r2 = C.evolve(A, B, 0.0, 0.0, pot, h, G=1e-10)
for tag, r in (("explicit (throat_evolution)", r1), ("implicit (condensate)", r2)):
    print(tag, "max|dphi|/v on diagonal:", np.nanmax(np.abs([r["phi"][k, k] - ph(k*h) for k in range(N+1)]))/v,
          " at t=0.25,1,4:", np.round([(r["phi"][int(t/h), int(t/h)] - ph(t))/v for t in (0.25, 1.0, 4.0)], 4))
