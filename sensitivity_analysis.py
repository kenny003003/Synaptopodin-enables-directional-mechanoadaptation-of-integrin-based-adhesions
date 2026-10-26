#!/usr/bin/env python3
"""
sensitivity_analysis.py
=======================
Parameter sensitivity of the orthogonal force-balance catch-slip detachment
model (Qu et al., "Synaptopodin enables directional mechanoadaptation of
integrin-based adhesions"). Reproduces Supplementary Tables 2-3 and the data
for Supplementary Figure 15.

Model (physical parameterization; see Supplementary Table 1 and Methods
Eqs. 7-16): N stress fibers around an elliptical cell apply force to integrin
bonds; each timestep a vertical displacement delta is solved so the bonded
fibers balance the applied perpendicular force. Each bond's normal (N_i) and
shear (S_i) load set a catch-slip off-rate
    k_off = k_off0 * exp[ (N_i/N0) * (Ks^n/(Ks^n+S_i^n)) * (1 + C*S_i^m/(K0^m+S_i^m)) ].
Genotypes differ only in C (WT 1.5, synaptopodin-deficient 2.7).

Outputs
    results_local_OAT.csv : local one-at-a-time elasticities S = (dY/Y)/(dtheta/theta),
                            uniform +/-10% perturbation, for WT and KO.
    results_sobol.csv     : variance-based Sobol first-order (S1) and total-effect (ST)
                            indices for WT detachment (Saltelli/Jansen, +/-30% box).

Run:  python3 sensitivity_analysis.py      (~5-8 min on a normal desktop)
"""
import numpy as np
import pandas as pd

RNG = np.random.default_rng(7)

AR = 1.4                                   # cell aspect ratio
PHYS = dict(a=1.0/np.sqrt(AR), b=np.sqrt(AR), r0=0.3, h=0.1,
            k0=1.0e6, F0=2.5, k_off0=0.8, N0=2.5, K_s=2.5, K_0=2.5,
            n_exp=1.1, m_exp=4, N=100, dt=0.01, nt=250, DET=0.05)
C_WT, C_KO = 1.5, 2.7


# --------------------------- vectorized model --------------------------------
def geom(P):
    N = int(P["N"])
    theta = np.linspace(0, 2*np.pi*(N-1)/N, N)
    r = P["a"]/np.sqrt(np.cos(theta)**2 + (P["a"]/P["b"])**2*np.sin(theta)**2)
    dr = r - P["r0"]
    L0 = np.sqrt(dr**2 + P["h"]**2)
    k_i = P["k0"]*P["h"]/L0
    return dr, L0, k_i


def solve_delta(delta, dr, L0, k_i, P, bonded, F, iters):
    """Damped-Newton force balance, vectorized over replicates. Uses the exact
    identities sin(phi)=hpd/L_f, cos(phi)=dr/L_f (phi = atan(hpd/dr))."""
    h, F0 = P["h"], P["F0"]
    for _ in range(iters):
        hpd = (h + delta)[:, None]
        L_f = np.sqrt(dr**2 + hpd**2)
        sinphi, cosphi = hpd/L_f, dr/L_f
        F_i = F0 + k_i*h*(L_f/L0 - 1.0)
        resid = np.sum(F_i*sinphi*bonded, axis=1) - F
        dF = k_i*h/L0*(hpd/L_f)
        dphi = dr/(dr**2 + hpd**2)
        dRv = np.sum((dF*sinphi + F_i*cosphi*dphi)*bonded, axis=1)
        dRv = np.where(np.abs(dRv) < 1e-12, 1e-12, dRv)
        delta = np.clip(delta - 0.3*(resid/dRv), 0.0, 20.0)
    return delta


def simulate(C, F, P, R, rng, first=35, warm=4):
    """Return the fraction of R Monte-Carlo cells that detach at applied force F."""
    N, nt, dt, DET = int(P["N"]), int(P["nt"]), P["dt"], P["DET"]
    h, F0, N0, k_off0 = P["h"], P["F0"], P["N0"], P["k_off0"]
    dr, L0, k_i = geom(P)
    Ksn, K0m = P["K_s"]**P["n_exp"], P["K_0"]**P["m_exp"]
    bonded = np.ones((R, N), dtype=bool)
    delta = np.zeros(R)
    detached = np.zeros(R, dtype=bool)
    stall = 0
    for t in range(nt):
        delta = solve_delta(delta, dr, L0, k_i, P, bonded, F, first if t == 0 else warm)
        hpd = (h + delta)[:, None]
        L_f = np.sqrt(dr**2 + hpd**2)
        sinphi, cosphi = hpd/L_f, dr/L_f
        F_i = F0 + k_i*h*(L_f/L0 - 1.0)
        N_i, S_i = F_i*sinphi, F_i*cosphi
        arg = (N_i/N0)*(Ksn/(Ksn + S_i**P["n_exp"]))*(1 + C*S_i**P["m_exp"]/(K0m + S_i**P["m_exp"]))
        koff = k_off0*np.exp(np.minimum(arg, 700.0))
        prob = 1 - np.exp(-koff*dt)
        broke = bonded & (prob >= rng.random((R, N)))
        bonded &= ~broke
        detached |= (bonded.mean(1) < DET)
        if detached.all():
            break
        stall = stall + 1 if not broke.any() else 0
        if stall >= 5:
            break
    return detached.mean()


# --------------------------- operating point ---------------------------------
def force_for_KO(target=0.50, lo=5.0, hi=320.0, reps=1500):
    for _ in range(14):
        mid = np.sqrt(lo*hi)
        if simulate(C_KO, mid, PHYS, reps, RNG) < target:
            lo = mid
        else:
            hi = mid
    return np.sqrt(lo*hi)

print("=== Operating point (C_WT=1.5, C_KO=2.7; force set so KO ~ 50% detachment) ===")
F_OP = force_for_KO()
wt0 = simulate(C_WT, F_OP, PHYS, 4000, RNG)
ko0 = simulate(C_KO, F_OP, PHYS, 4000, RNG)
print(f"  F_OP = {F_OP:.2f}   WT = {100*wt0:.1f}%   KO = {100*ko0:.1f}%   gap = {100*(ko0-wt0):+.1f} pp\n")


# --------------------------- local one-at-a-time -----------------------------
print("=== Local one-at-a-time sensitivity (uniform +/-10%) ===")
R = 4000
params = ["a", "b", "r0", "h", "k0", "F0", "k_off0", "N0", "K_s", "K_0", "n_exp", "m_exp", "N"]
rows = []
for p in params:
    th = PHYS[p]
    loP, hiP = dict(PHYS), dict(PHYS)
    loP[p], hiP[p] = th*0.9, th*1.1
    if p == "N":                                   # keep bond count an integer
        loP[p], hiP[p] = int(round(th*0.9)), int(round(th*1.1))
    dth = hiP[p] - loP[p]
    wL, wH = simulate(C_WT, F_OP, loP, R, RNG), simulate(C_WT, F_OP, hiP, R, RNG)
    kL, kH = simulate(C_KO, F_OP, loP, R, RNG), simulate(C_KO, F_OP, hiP, R, RNG)
    S_WT = ((wH - wL)/wt0)/(dth/th) if wt0 else 0.0
    S_KO = ((kH - kL)/ko0)/(dth/th) if ko0 else 0.0
    rows.append(dict(param=p, S_WT=S_WT, S_KO=S_KO))
dfO = pd.DataFrame(rows).sort_values("S_WT", key=lambda s: s.abs(), ascending=False)
dfO.to_csv("results_local_OAT.csv", index=False)
print(dfO.to_string(index=False, float_format=lambda v: f"{v:+.2f}"))
print("  -> results_local_OAT.csv\n")


# --------------------------- Sobol indices -----------------------------------
print("=== Variance-based Sobol indices (Saltelli/Jansen, +/-30% box, WT) ===")
sob = ["a", "b", "r0", "h", "k0", "F0", "k_off0", "N0", "K_s", "K_0", "n_exp", "m_exp"]
d = len(sob)
lo = np.array([PHYS[p]*0.7 for p in sob])
hi = np.array([PHYS[p]*1.3 for p in sob])
Nb, REP = 48, 200

def model_vec(row):
    P = dict(PHYS)
    for j, p in enumerate(sob):
        P[p] = row[j]
    return simulate(C_WT, F_OP, P, REP, RNG)

A = lo + RNG.random((Nb, d))*(hi - lo)
B = lo + RNG.random((Nb, d))*(hi - lo)
fA = np.array([model_vec(A[i]) for i in range(Nb)])
fB = np.array([model_vec(B[i]) for i in range(Nb)])
VarY = np.var(np.concatenate([fA, fB]), ddof=1)
S1 = np.zeros(d); ST = np.zeros(d)
for j in range(d):
    AB = A.copy(); AB[:, j] = B[:, j]
    fAB = np.array([model_vec(AB[i]) for i in range(Nb)])
    S1[j] = np.mean(fB*(fAB - fA))/VarY if VarY > 0 else 0.0
    ST[j] = np.mean((fA - fAB)**2)/(2*VarY) if VarY > 0 else 0.0
dfS = pd.DataFrame(dict(param=sob, S1=S1, ST=ST)).sort_values("ST", ascending=False)
dfS.to_csv("results_sobol.csv", index=False)
print(dfS.to_string(index=False, float_format=lambda v: f"{v:+.3f}"))
print("  -> results_sobol.csv")
print(f"\nDONE. Operating point F={F_OP:.1f}, WT {100*wt0:.0f}% / KO {100*ko0:.0f}%.")
