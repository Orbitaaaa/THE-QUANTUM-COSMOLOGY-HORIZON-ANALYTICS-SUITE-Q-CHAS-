import numpy as np, time, warnings
warnings.simplefilter("ignore")
from babyuniverse import condensate_eft as E
from scipy.integrate import solve_ivp
# reference: per-mode solve_ivp (the old implementation) on a small case
def ref(mu, alpha, eps0, J, sigma, L, N, t_end, seed_amp=1e-8, seed=0):
    rng = np.random.default_rng(seed)
    k = 2 * np.pi * np.fft.rfftfreq(N, d=L / N)
    phase = rng.uniform(0, 2 * np.pi, len(k))
    amp = np.where((k > 0) & (k < k.max() / 2), 1.0, 0.0)
    pi0 = np.fft.irfft(amp * np.exp(1j * phase), n=N)
    c0 = np.fft.rfft(pi0 * seed_amp / max(np.std(pi0), 1e-300))
    C = alpha / mu**2; A = 2 + eps0
    eps_t = (lambda t: eps0 * np.exp(-t)) if sigma == 0 else (lambda t: -sigma + (eps0 + sigma) * np.exp(-t))
    ts = np.linspace(0, t_end, 51); pw = np.zeros(51)
    for kk, ck in zip(k, c0):
        if abs(ck) == 0: continue
        f = lambda t, z: [z[1], -z[1] - (eps_t(t) * (kk*np.exp(-t))**2 + 0.5*C*(kk*np.exp(-t))**4 - J*(kk*np.exp(-t))**2)/A*z[0]]
        sol = solve_ivp(f, (0, t_end), [1.0, 0.0], t_eval=ts, rtol=1e-10, atol=1e-13, max_step=0.01)
        pw += 2*abs(ck)**2*sol.y[0]**2
    rms = np.sqrt(pw)/N
    return rms.max()/rms[0]
for (mu, al, e0, J, sg, L, t_end) in [(30.0, 1.0, -0.05, 0.0, 0.0, 6.0, 9.6), (30.0, 1.0, -0.05, 0.0, 0.05, 6.0, 20.0), (50.0, 1.0, -0.02, 0.01, 0.0, 10.0, 9.6)]:
    t0 = time.time(); a = E.linear_rms_1p1(mu, al, e0, J=J, sigma=sg, L=L, N=128, t_end=t_end)["amplification"]; t1 = time.time()
    b = ref(mu, al, e0, J, sg, L, 128, t_end); t2 = time.time()
    print(f"new {a:.6f} ({t1-t0:.2f}s)  ref {b:.6f} ({t2-t1:.2f}s)  rel {abs(a/b-1):.1e}")
